from pathlib import Path
import click
import numpy as np
from typing import List, Tuple, Dict, Optional
import json
from collections import defaultdict
import random
from dataclasses import dataclass
import sys
from loguru import logger

@dataclass
class SequencePattern:
    """Store information about detected sequence patterns."""
    start: int
    length: int
    pattern_type: str  # 'UMI' or 'BARCODE'
    entropy: float
    coverage: int
    base_frequencies: Dict[str, Dict[str, float]]  # Position -> base -> frequency
    consensus: str = ""  # For barcodes only

class FastqReader:
    """Memory-efficient FASTQ reader with reservoir sampling."""
    
    @staticmethod
    def reservoir_sample(fastq_path: str, n: int = 1000, 
                        seed: Optional[int] = None) -> List[str]:
        """
        Reservoir sampling of FASTQ reads.
        
        Args:
            fastq_path: Path to FASTQ file
            n: Number of reads to sample
            seed: Random seed for reproducibility
        """
        if seed is not None:
            random.seed(seed)
            
        reservoir = []
        
        try:
            with open(fastq_path, 'r') as f:
                # Process first n reads
                for _ in range(n):
                    try:
                        # Read 4 lines (one FASTQ record)
                        lines = [next(f) for _ in range(4)]
                        reservoir.append(lines[1].strip())  # Sequence line
                    except StopIteration:
                        return reservoir
                    
                # Reservoir sampling for remaining reads
                for i, line in enumerate(f, start=n*4):
                    if i % 4 == 1:  # Only process sequence lines
                        j = random.randint(0, i//4)
                        if j < n:
                            reservoir[j] = line.strip()
        except Exception as e:
            logger.error(f"Error reading FASTQ file: {e}")
            sys.exit(1)
            
        return reservoir

def calculate_base_frequencies(bases: List[str]) -> Dict[str, float]:
    """Calculate frequency of each nucleotide."""
    freqs = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    total = len(bases)
    
    for base in bases:
        if base in freqs:
            freqs[base] += 1
    
    # Normalize frequencies
    return {base: count/total for base, count in freqs.items()}

def calculate_entropy(frequencies: Dict[str, float]) -> float:
    """Calculate Shannon entropy from base frequencies."""
    return -sum(p * np.log2(p) if p > 0 else 0 for p in frequencies.values())

def calculate_position_stats(sequences: List[str], 
                           min_coverage: int = 100,
                           from_end: bool = False) -> Dict:
    """
    Calculate positional statistics from either 5' or 3' end.
    
    Args:
        sequences: List of sequence strings
        min_coverage: Minimum number of reads required for analysis
        from_end: If True, calculate from 3' end
        
    Returns:
        Dict containing position-wise statistics
    """
    stats = defaultdict(dict)
    
    if from_end:
        sequences = [seq[::-1] for seq in sequences]  # Reverse sequences for 3' analysis
    
    max_len = max(len(seq) for seq in sequences)
    
    for pos in range(max_len):
        bases = [seq[pos] for seq in sequences if len(seq) > pos]
        if len(bases) < min_coverage:
            break
            
        freqs = calculate_base_frequencies(bases)
        entropy = calculate_entropy(freqs)
        
        stats[pos] = {
            'frequencies': freqs,
            'entropy': entropy,
            'coverage': len(bases)
        }
    
    return stats

def detect_patterns(stats: Dict, 
                   window_size: int = 6,
                   umi_entropy_threshold: float = 1.9,
                   barcode_entropy_threshold: float = 1.5) -> List[SequencePattern]:
    """
    Detect UMIs and barcodes based on positional statistics.
    
    Args:
        stats: Position-wise statistics
        window_size: Size of window for pattern detection
        umi_entropy_threshold: Minimum entropy for UMI detection
        barcode_entropy_threshold: Maximum entropy for barcode detection
    """
    patterns = []
    positions = sorted(stats.keys())
    
    if len(positions) < window_size:
        return patterns
    
    for start_pos in positions[:-window_size + 1]:
        window_stats = [stats[pos]['entropy'] 
                       for pos in range(start_pos, start_pos + window_size)]
        window_freqs = {
            pos: stats[pos]['frequencies']
            for pos in range(start_pos, start_pos + window_size)
        }
        
        avg_entropy = np.mean(window_stats)
        min_coverage = min(stats[pos]['coverage'] 
                         for pos in range(start_pos, start_pos + window_size))
        
        if avg_entropy > umi_entropy_threshold:
            patterns.append(SequencePattern(
                start=start_pos,
                length=window_size,
                pattern_type='UMI',
                entropy=avg_entropy,
                coverage=min_coverage,
                base_frequencies=window_freqs
            ))
        elif avg_entropy < barcode_entropy_threshold:
            # Generate consensus for barcode
            consensus = ""
            for pos in range(start_pos, start_pos + window_size):
                freqs = stats[pos]['frequencies']
                max_base = max(freqs.items(), key=lambda x: x[1])[0]
                consensus += max_base
                
            patterns.append(SequencePattern(
                start=start_pos,
                length=window_size,
                pattern_type='BARCODE',
                entropy=avg_entropy,
                coverage=min_coverage,
                base_frequencies=window_freqs,
                consensus=consensus
            ))
            
    return patterns

def generate_umitools_config(five_prime_patterns: List[SequencePattern],
                           three_prime_patterns: List[SequencePattern],
                           output_path: str):
    """Generate UMI-tools compatible configuration."""
    umi_patterns = {
        'five_prime': [p for p in five_prime_patterns if p.pattern_type == 'UMI'],
        'three_prime': [p for p in three_prime_patterns if p.pattern_type == 'UMI']
    }
    
    if not any(umi_patterns.values()):
        logger.info("No UMI patterns detected - skipping UMI-tools config generation")
        return
    
    config = {
        'extract_method': 'string',
        'pattern': ''
    }
    
    # Build pattern string based on detected patterns
    pattern_parts = []
    
    # Add 5' UMI if present
    if umi_patterns['five_prime']:
        p = umi_patterns['five_prime'][0]
        if p.start > 0:
            pattern_parts.append(f".{{{p.start}}}")
        pattern_parts.append(f"(?P<umi_1>.{{{p.length}}})")
        logger.debug(f"Added 5' UMI pattern of length {p.length}")
    
    # Add middle part if we have both 5' and 3' UMIs
    if umi_patterns['five_prime'] and umi_patterns['three_prime']:
        pattern_parts.append(".*")
    
    # Add 3' UMI if present
    if umi_patterns['three_prime']:
        p = umi_patterns['three_prime'][0]
        pattern_parts.append(f"(?P<umi_2>.{{{p.length}}})")
        logger.debug(f"Added 3' UMI pattern of length {p.length}")
    
    config['pattern'] = ''.join(pattern_parts)
    
    logger.info(f"Generated UMI-tools pattern: {config['pattern']}")
    
    with open(output_path, 'w') as f:
        json.dump(config, f, indent=2)

@click.command()
@click.argument('fastq_path', type=click.Path(exists=True))
@click.option('--sample-size', '-n', default=10000,
              help='Number of reads to sample')
@click.option('--window-size', '-w', default=6,
              help='Window size for pattern detection')
@click.option('--min-coverage', '-c', default=100,
              help='Minimum read coverage for analysis')
@click.option('--output-prefix', '-o', default='pattern_analysis',
              help='Prefix for output files')
@click.option('--seed', '-s', default=None, type=int,
              help='Random seed for reproducibility')
def main(fastq_path: str, sample_size: int, window_size: int, 
         min_coverage: int, output_prefix: str, seed: Optional[int]):
    """
    Analyze FASTQ file for UMIs and barcodes.
    
    This tool performs statistical analysis of read patterns to identify potential
    UMIs and barcodes. It analyzes both 5' and 3' ends of reads using entropy-based
    detection methods.
    """
    logger.info(f"Starting analysis of {fastq_path}")
    
    # Read sample of sequences
    reader = FastqReader()
    sequences = reader.reservoir_sample(fastq_path, n=sample_size, seed=seed)
    logger.info(f"Sampled {len(sequences)} sequences")
    
    # Calculate statistics
    logger.info("Calculating position-wise statistics")
    five_prime_stats = calculate_position_stats(
        sequences, min_coverage=min_coverage, from_end=False
    )
    three_prime_stats = calculate_position_stats(
        sequences, min_coverage=min_coverage, from_end=True
    )
    
    # Detect patterns
    logger.info("Detecting sequence patterns")
    five_prime_patterns = detect_patterns(
        five_prime_stats, window_size=window_size
    )
    three_prime_patterns = detect_patterns(
        three_prime_stats, window_size=window_size
    )
    
    # Generate outputs
    logger.info("Generating output files")
    results = {
        'parameters': {
            'sample_size': sample_size,
            'window_size': window_size,
            'min_coverage': min_coverage,
            'seed': seed
        },
        'five_prime_patterns': [vars(p) for p in five_prime_patterns],
        'three_prime_patterns': [vars(p) for p in three_prime_patterns]
    }
    
    # Write results
    output_path = Path(output_prefix)
    with open(f"{output_path}_results.json", 'w') as f:
        json.dump(results, f, indent=2)
    logger.info(f"Written detailed results to {output_path}_results.json")
    
    # Generate UMI-tools compatible output if UMIs detected
    if any(p.pattern_type == 'UMI' for p in five_prime_patterns + three_prime_patterns):
        generate_umitools_config(
            five_prime_patterns, 
            three_prime_patterns,
            f"{output_path}_umitools_config.json"
        )
        logger.info(f"Written UMI-tools config to {output_path}_umitools_config.json")
    
    logger.info("Analysis complete")

if __name__ == '__main__':
    main()