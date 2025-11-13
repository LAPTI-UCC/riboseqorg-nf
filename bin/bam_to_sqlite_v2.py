#!/usr/bin/env python3
"""
Convert a BAM file to a sqlite read file for ribosome profiling analysis.

This refactored version improves on the original by:
- Using classes for better organization
- Reducing redundant lookups with caching
- Improving readability with better naming
- Adding comprehensive logging
- Optimizing data structure operations
"""

import sys
import logging
from typing import Dict, Set, Tuple, Optional, List
from dataclasses import dataclass, field
from collections import defaultdict, Counter
import argparse

import pysam
import sqlite3
from sqlitedict import SqliteDict


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class TranscriptInfo:
    """Stores information about a transcript."""
    chrom: str
    strand: str
    cds_start: Optional[int]
    cds_stop: Optional[int]
    length: int
    tran_type: int
    exons: List[Tuple[int, int]] = field(default_factory=list)


@dataclass
class ReadStats:
    """Aggregates read statistics."""
    total_reads: int = 0
    mapped_reads: int = 0
    unmapped_reads: int = 0
    ambiguous_reads: int = 0
    unambiguous_coding: int = 0
    unambiguous_noncoding: int = 0

    # Read length distributions
    read_lengths: Counter = field(default_factory=Counter)
    unambig_read_lengths: Counter = field(default_factory=Counter)

    # Nucleotide composition
    nuc_composition_mapped: Dict = field(default_factory=lambda: defaultdict(lambda: defaultdict(Counter)))
    nuc_composition_unmapped: Dict = field(default_factory=lambda: defaultdict(lambda: defaultdict(Counter)))
    threeprime_nuc_mapped: Dict = field(default_factory=lambda: defaultdict(lambda: defaultdict(Counter)))
    threeprime_nuc_unmapped: Dict = field(default_factory=lambda: defaultdict(lambda: defaultdict(Counter)))
    dinuc_counts: Dict = field(default_factory=lambda: defaultdict(Counter))

    # Unmapped read tracking
    unmapped_seqs: Counter = field(default_factory=Counter)


class TranscriptomeDB:
    """Handles loading and querying transcriptome annotation data."""

    def __init__(self, db_path: str):
        self.transcripts: Dict[str, TranscriptInfo] = {}
        self._load_from_db(db_path)

    def _load_from_db(self, db_path: str):
        """Load transcript information from SQLite database."""
        logger.info(f"Loading transcriptome annotation from {db_path}")

        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Load transcript information
        cursor.execute("SELECT transcript, cds_start, cds_stop, length, strand, chrom, tran_type FROM transcripts")
        for row in cursor.fetchall():
            tran_id = str(row[0])
            self.transcripts[tran_id] = TranscriptInfo(
                chrom=row[5],
                strand=row[4],
                cds_start=row[1],
                cds_stop=row[2],
                length=row[3],
                tran_type=row[6],
                exons=[]
            )

        # Load exon information
        cursor.execute("SELECT * FROM exons")
        for row in cursor.fetchall():
            tran_id = str(row[0])
            if tran_id in self.transcripts:
                self.transcripts[tran_id].exons.append((row[1], row[2]))

        conn.close()
        logger.info(f"Loaded {len(self.transcripts)} transcripts")

    def get(self, transcript_id: str) -> Optional[TranscriptInfo]:
        """Get transcript info by ID."""
        return self.transcripts.get(transcript_id)

    def transcript_to_genomic_pos(self, transcript_id: str, pos: int) -> Optional[Tuple[str, int]]:
        """Convert transcript position to genomic coordinates."""
        tran_info = self.get(transcript_id)
        if not tran_info:
            return None

        exons = sorted(tran_info.exons)

        if tran_info.strand == "+":
            for exon_start, exon_end in exons:
                exon_len = exon_end - exon_start
                if pos > exon_len:
                    pos = (pos - exon_len) - 1
                else:
                    genomic_pos = (exon_start + pos) - 1
                    return (tran_info.chrom, genomic_pos)
        else:  # strand == "-"
            for exon_start, exon_end in reversed(exons):
                exon_len = exon_end - exon_start
                if pos > exon_len:
                    pos = (pos - exon_len) - 1
                else:
                    genomic_pos = (exon_end - pos) + 1
                    return (tran_info.chrom, genomic_pos)

        return None


class ReadProcessor:
    """Processes ribosome profiling reads and generates statistics."""

    def __init__(self, transcriptome_db: TranscriptomeDB):
        self.transcriptome_db = transcriptome_db
        self.stats = ReadStats()
        self.missing_transcripts: Set[str] = set()

        # Per-transcript read data
        self.transcript_reads: Dict[str, Dict] = defaultdict(
            lambda: {
                "ambig": defaultdict(lambda: defaultdict(int)),
                "unambig": defaultdict(lambda: defaultdict(int)),
                "mismatches": {},
                "seq": defaultdict(lambda: defaultdict(int))
            }
        )

        # Metagene analysis data
        self.trip_periodicity = {
            "fiveprime": defaultdict(lambda: {"0": 0.0, "1": 0.0, "2": 0.0}),
            "threeprime": defaultdict(lambda: {"0": 0.0, "1": 0.0, "2": 0.0})
        }
        self.metagene_start = {
            "fiveprime": defaultdict(lambda: defaultdict(int)),
            "threeprime": defaultdict(lambda: defaultdict(int))
        }
        self.metagene_stop = {
            "fiveprime": defaultdict(lambda: defaultdict(int)),
            "threeprime": defaultdict(lambda: defaultdict(int))
        }

    @staticmethod
    def get_read_count(read_name: str) -> int:
        """Extract read count from collapsed read names (e.g., read100_x100)."""
        try:
            return int(read_name.split('_x')[-1])
        except (IndexError, ValueError):
            return 1

    @staticmethod
    def normalize_transcript_id(tran_id: str) -> str:
        """Normalize transcript ID by removing special characters."""
        return tran_id.replace("-", "_").replace("(", "").replace(")", "")

    def parse_mismatches(self, md_tag: str, read_seq: str) -> Dict[int, str]:
        """Parse MD tag to extract mismatch positions and bases."""
        mismatches = {}
        nucs = set("ATGC")
        total_so_far = 0
        prev_digits = ""

        for char in md_tag:
            if char in nucs:
                if prev_digits:
                    total_so_far += int(prev_digits)
                    prev_digits = ""
                pos = total_so_far + len(mismatches)
                if pos < len(read_seq):
                    mismatches[pos] = read_seq[pos]
            elif char not in ("^", "N"):  # Skip deletions and N's
                prev_digits += char

        return mismatches

    def analyze_nucleotide_composition(self, seq: str, read_len: int,
                                      read_count: int, is_mapped: bool):
        """Analyze nucleotide composition of reads."""
        nuc_dict = self.stats.nuc_composition_mapped if is_mapped else self.stats.nuc_composition_unmapped
        threeprime_dict = self.stats.threeprime_nuc_mapped if is_mapped else self.stats.threeprime_nuc_unmapped

        # 5' to 3' nucleotide composition
        for i, base in enumerate(seq):
            nuc_dict[read_len][i][base] += read_count

        # 3' end nucleotide composition
        for i in range(len(seq), 0, -1):
            dist = i - len(seq)
            threeprime_dict[read_len][dist][seq[dist]] += read_count

        # Dinucleotide composition (only for mapped reads)
        if is_mapped:
            for i in range(len(seq) - 1):
                dinuc = seq[i:i+2]
                if len(dinuc) == 2:
                    self.stats.dinuc_counts[read_len][dinuc] += read_count

    def process_read_alignments(self, read_name: str, alignments: List[Tuple[str, int, List]]):
        """
        Process all alignments for a single read.

        Args:
            read_name: Name of the read
            alignments: List of (transcript_id, position, tags) tuples
        """
        if not alignments:
            return

        read_count = self.get_read_count(read_name)
        read_seq = alignments[0][3] if len(alignments[0]) > 3 else ""
        read_len = len(read_seq)

        # Get unique genomic positions
        genomic_positions = []
        for tran_id, pos, _ in alignments[:3]:  # Use first 3 items
            tran_id = self.normalize_transcript_id(tran_id)
            genomic_pos = self.transcriptome_db.transcript_to_genomic_pos(tran_id, pos)

            if genomic_pos is None:
                self.missing_transcripts.add(tran_id)
                continue

            if genomic_pos not in genomic_positions:
                genomic_positions.append(genomic_pos)

        if not genomic_positions:
            return

        is_unambiguous = len(genomic_positions) == 1

        if is_unambiguous:
            self.stats.unambig_read_lengths[read_len] += read_count

        # Process each alignment
        is_coding = False
        for tran_id, pos, tags, *_ in alignments:
            tran_id = self.normalize_transcript_id(tran_id)
            tran_info = self.transcriptome_db.get(tran_id)

            if not tran_info:
                self.missing_transcripts.add(tran_id)
                continue

            # Store read position
            if is_unambiguous:
                self.transcript_reads[tran_id]["unambig"][read_len][pos] += read_count

                # Check if read is in coding region
                if tran_info.cds_start and tran_info.cds_stop:
                    if tran_info.cds_start <= pos <= tran_info.cds_stop:
                        is_coding = True
            else:
                self.transcript_reads[tran_id]["ambig"][read_len][pos] += read_count

            # Process mismatches
            nm_tag = dict(tags).get("NM", 0)
            if nm_tag > 0:
                md_tag = dict(tags).get("MD", "")
                if md_tag:
                    mismatches = self.parse_mismatches(md_tag, read_seq)
                    for mismatch_offset, base in mismatches.items():
                        mismatch_pos = pos + mismatch_offset
                        self.transcript_reads[tran_id]["seq"][mismatch_pos][base] += read_count

        # Update coding/noncoding stats
        if is_unambiguous:
            if is_coding:
                self.stats.unambiguous_coding += read_count
            else:
                self.stats.unambiguous_noncoding += read_count
        else:
            self.stats.ambiguous_reads += read_count

    def compute_metagene_analysis(self):
        """Compute triplet periodicity and metagene profiles."""
        logger.info("Computing metagene analysis and triplet periodicity")

        for tran_id, read_data in self.transcript_reads.items():
            tran_info = self.transcriptome_db.get(tran_id)
            if not tran_info:
                continue

            cds_start = tran_info.cds_start
            cds_stop = tran_info.cds_stop

            # Skip transcripts without proper CDS or UTRs
            if not cds_start or not cds_stop:
                continue
            if cds_start <= 1 or cds_stop >= tran_info.length:
                continue
            if tran_info.tran_type != 1:
                continue

            # Process unambiguous reads
            for read_len, positions in read_data["unambig"].items():
                for pos, count in positions.items():
                    for prime_type in ["fiveprime", "threeprime"]:
                        if prime_type == "fiveprime":
                            real_pos = pos - cds_start
                            rel_stop_pos = pos - cds_stop
                        else:
                            real_pos = (pos + read_len) - cds_start
                            rel_stop_pos = (pos + read_len) - cds_stop

                        # Triplet periodicity (within CDS)
                        if cds_start <= pos <= cds_stop:
                            frame = str(real_pos % 3)
                            self.trip_periodicity[prime_type][read_len][frame] += count

                        # Metagene profile around start codon
                        if -600 < real_pos < 601:
                            self.metagene_start[prime_type][read_len][real_pos] += count

                        # Metagene profile around stop codon
                        if -600 < rel_stop_pos < 601:
                            self.metagene_stop[prime_type][read_len][rel_stop_pos] += count

    def compute_offsets(self) -> Dict:
        """Compute optimal offsets for each read length."""
        logger.info("Computing P-site offsets")

        offsets = {
            "fiveprime": {"offsets": {}, "read_scores": {}},
            "threeprime": {"offsets": {}, "read_scores": {}}
        }

        for prime_type in ["fiveprime", "threeprime"]:
            # Compute triplet periodicity scores
            for read_len, frames in self.trip_periodicity[prime_type].items():
                counts = sorted([frames["0"], frames["1"], frames["2"]], reverse=True)
                if counts[0] > 0:
                    score = counts[1] / counts[0] if counts[0] > 0 else 0.0
                    offsets[prime_type]["read_scores"][read_len] = score
                else:
                    offsets[prime_type]["read_scores"][read_len] = 0.0

            # Compute offsets from metagene profiles
            for read_len, positions in self.metagene_start[prime_type].items():
                max_pos = 0
                max_count = 0

                for pos, count in positions.items():
                    # Only consider positions within reasonable range
                    if 10 <= abs(pos) <= (read_len - 10):
                        if count > max_count:
                            max_count = count
                            max_pos = pos

                if prime_type == "fiveprime":
                    # -3 to get from P-site to A-site, +1 for 1-based coords = -2
                    offsets[prime_type]["offsets"][read_len] = abs(max_pos - 2)
                else:
                    # +3 to get from P-site to A-site, -1 for 1-based coords = +2
                    offsets[prime_type]["offsets"][read_len] = (max_pos * -1) + 2

        return offsets

    def compute_region_totals(self, offsets: Dict) -> Dict:
        """Compute read totals for 5'UTR, CDS, and 3'UTR regions."""
        logger.info("Computing region-specific read totals")

        totals = {
            "unambiguous_all": {},
            "unambiguous_fiveprime": {},
            "unambiguous_cds": {},
            "unambiguous_threeprime": {},
            "ambiguous_all": {},
            "ambiguous_fiveprime": {},
            "ambiguous_cds": {},
            "ambiguous_threeprime": {}
        }

        for tran_id, read_data in self.transcript_reads.items():
            tran_info = self.transcriptome_db.get(tran_id)
            if not tran_info:
                continue

            five_total = cds_total = three_total = 0
            ambig_five = ambig_cds = ambig_three = 0

            cds_start = tran_info.cds_start
            cds_stop = tran_info.cds_stop

            # Process unambiguous reads
            for read_len, positions in read_data["unambig"].items():
                offset = offsets["fiveprime"]["offsets"].get(read_len, 15)

                for pos, count in positions.items():
                    real_pos = pos + offset

                    if not cds_start or not cds_stop:
                        three_total += count
                    elif real_pos < cds_start:
                        five_total += count
                    elif cds_start <= real_pos <= cds_stop:
                        cds_total += count
                    else:
                        three_total += count

            # Process ambiguous reads
            for read_len, positions in read_data["ambig"].items():
                offset = offsets["fiveprime"]["offsets"].get(read_len, 15)

                for pos, count in positions.items():
                    real_pos = pos + offset

                    if not cds_start or not cds_stop:
                        ambig_three += count
                    elif real_pos < cds_start:
                        ambig_five += count
                    elif cds_start <= real_pos <= cds_stop:
                        ambig_cds += count
                    else:
                        ambig_three += count

            # Store totals
            totals["unambiguous_all"][tran_id] = five_total + cds_total + three_total
            totals["unambiguous_fiveprime"][tran_id] = five_total
            totals["unambiguous_cds"][tran_id] = cds_total
            totals["unambiguous_threeprime"][tran_id] = three_total

            totals["ambiguous_all"][tran_id] = five_total + cds_total + three_total + ambig_five + ambig_cds + ambig_three
            totals["ambiguous_fiveprime"][tran_id] = five_total + ambig_five
            totals["ambiguous_cds"][tran_id] = cds_total + ambig_cds
            totals["ambiguous_threeprime"][tran_id] = three_total + ambig_three

        return totals


def process_bam_file(bam_path: str, annotation_db_path: str, output_path: str):
    """Main function to process BAM file and generate SQLite output."""
    logger.info(f"Processing BAM file: {bam_path}")

    # Load transcriptome annotation
    transcriptome_db = TranscriptomeDB(annotation_db_path)

    # Initialize processor
    processor = ReadProcessor(transcriptome_db)

    # Open BAM file
    pysam.set_verbosity(0)
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Check if BAM is sorted by read name
    header = bam.header.get("HD", {})
    if header.get("SO") != "queryname":
        logger.error("BAM file must be sorted by read name. Use: samtools sort -n input.bam -o output.bam")
        sys.exit(1)

    # Process reads
    logger.info("Processing reads from BAM file")
    current_read_name = None
    current_alignments = []
    reads_processed = 0

    for read in bam.fetch(until_eof=True):
        read_count = processor.get_read_count(read.query_name)
        processor.stats.total_reads += read_count

        if read.is_unmapped:
            processor.stats.unmapped_reads += read_count
            seq = read.query_sequence
            processor.stats.unmapped_seqs[seq] += read_count
            processor.analyze_nucleotide_composition(seq, len(seq), read_count, is_mapped=False)
            continue

        processor.stats.mapped_reads += read_count

        # Get alignment info
        tran_id = bam.get_reference_name(read.reference_id).split(".")[0]
        pos = read.reference_start
        seq = read.query_sequence
        tags = read.get_tags()

        # Group alignments by read name
        if read.query_name != current_read_name:
            if current_alignments:
                processor.process_read_alignments(current_read_name, current_alignments)
                reads_processed += 1
                if reads_processed % 100000 == 0:
                    logger.info(f"Processed {reads_processed:,} reads")

            current_read_name = read.query_name
            current_alignments = [(tran_id, pos, tags, seq)]

            # Analyze nucleotide composition once per read
            read_len = len(seq)
            processor.stats.read_lengths[read_len] += read_count
            processor.analyze_nucleotide_composition(seq, read_len, read_count, is_mapped=True)
        else:
            current_alignments.append((tran_id, pos, tags, seq))

    # Process last read
    if current_alignments:
        processor.process_read_alignments(current_read_name, current_alignments)

    bam.close()

    logger.info(f"Total reads: {processor.stats.total_reads:,}")
    logger.info(f"Mapped reads: {processor.stats.mapped_reads:,}")
    logger.info(f"Unmapped reads: {processor.stats.unmapped_reads:,}")
    logger.info(f"Ambiguous reads: {processor.stats.ambiguous_reads:,}")

    # Warn about missing transcripts
    if processor.missing_transcripts:
        logger.warning(f"{len(processor.missing_transcripts)} transcript(s) found in BAM but not in annotation")
        logger.warning(f"Examples: {sorted(list(processor.missing_transcripts))[:10]}")

    # Compute metagene analysis
    processor.compute_metagene_analysis()

    # Compute offsets
    offsets = processor.compute_offsets()

    # Compute region totals
    region_totals = processor.compute_region_totals(offsets)

    # Write output
    logger.info(f"Writing results to {output_path}")
    write_output_db(output_path, processor, offsets, region_totals)

    logger.info("Processing complete")


def write_output_db(output_path: str, processor: ReadProcessor,
                   offsets: Dict, region_totals: Dict):
    """Write results to SQLite database."""
    db = SqliteDict(output_path, autocommit=False)

    # Write transcript read data
    for tran_id, data in processor.transcript_reads.items():
        db[tran_id] = dict(data)

    # Write statistics
    db["mapped_reads"] = processor.stats.mapped_reads
    db["unmapped_reads"] = processor.stats.unmapped_reads
    db["offsets"] = offsets
    db["trip_periodicity"] = dict(processor.trip_periodicity)
    db["nuc_counts"] = {
        "mapped": {k: dict(v) for k, v in processor.stats.nuc_composition_mapped.items()},
        "unmapped": {k: dict(v) for k, v in processor.stats.nuc_composition_unmapped.items()}
    }
    db["dinuc_counts"] = {k: dict(v) for k, v in processor.stats.dinuc_counts.items()}
    db["threeprime_nuc_counts"] = {
        "mapped": {k: dict(v) for k, v in processor.stats.threeprime_nuc_mapped.items()},
        "unmapped": {k: dict(v) for k, v in processor.stats.threeprime_nuc_unmapped.items()}
    }
    db["metagene_counts"] = {
        "fiveprime": {k: dict(v) for k, v in processor.metagene_start["fiveprime"].items()},
        "threeprime": {k: dict(v) for k, v in processor.metagene_start["threeprime"].items()}
    }
    db["stop_metagene_counts"] = {
        "fiveprime": {k: dict(v) for k, v in processor.metagene_stop["fiveprime"].items()},
        "threeprime": {k: dict(v) for k, v in processor.metagene_stop["threeprime"].items()}
    }
    db["read_lengths"] = dict(processor.stats.read_lengths)
    db["unambig_read_lengths"] = dict(processor.stats.unambig_read_lengths)
    db["coding_counts"] = processor.stats.unambiguous_coding
    db["noncoding_counts"] = processor.stats.unambiguous_noncoding
    db["ambiguous_counts"] = processor.stats.ambiguous_reads
    db["frequent_unmapped_reads"] = processor.stats.unmapped_seqs.most_common(2000)
    db["description"] = "Null"
    db["cutadapt_removed"] = 0
    db["rrna_removed"] = 0
    db["removed_minus_m"] = 0

    # Write region totals
    for key, values in region_totals.items():
        db[key + "_totals"] = values

    db.commit()
    db.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Convert a BAM file to ribosome profiling SQLite format (refactored version).'
    )
    parser.add_argument('--bam', required=True, help='Path to name-sorted BAM file')
    parser.add_argument('--annotation', required=True, help='Path to annotation SQLite file')
    parser.add_argument('--output', required=True, help='Path to output SQLite file')

    args = parser.parse_args()

    process_bam_file(args.bam, args.annotation, args.output)
