#!/usr/bin/env python3
"""
Process seqspec YAML files to extract adapter and barcode information
for downstream processing in the riboseq pipeline.
"""

import argparse
import json
import sys
import yaml
from typing import Dict, List, Optional, Any


def parse_seqspec(seqspec_file: str) -> Dict[str, Any]:
    """Parse a seqspec YAML file."""
    with open(seqspec_file, 'r') as f:
        return yaml.safe_load(f)


def extract_region_info(region: Dict[str, Any], parent_info: Optional[Dict] = None) -> Dict[str, Any]:
    """
    Recursively extract information about a region in the seqspec.

    Args:
        region: Region dictionary from seqspec
        parent_info: Information from parent region (for position tracking)

    Returns:
        Dictionary with region information
    """
    info = {
        'region_id': region.get('region_id', ''),
        'region_type': region.get('region_type', ''),
        'name': region.get('name', ''),
        'sequence_type': region.get('sequence_type', ''),
        'sequence': region.get('sequence', ''),
        'min_len': region.get('min_len', 0),
        'max_len': region.get('max_len', 0),
        'onlist': region.get('onlist', None),
    }

    # Handle nested regions
    if 'regions' in region and region['regions']:
        info['subregions'] = []
        for subregion in region['regions']:
            info['subregions'].append(extract_region_info(subregion, info))

    return info


def find_adapters(seqspec: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Find all adapter sequences in the seqspec.

    Returns:
        List of adapter information dictionaries
    """
    adapters = []

    def search_regions(regions: List[Dict], position: str = ""):
        for i, region in enumerate(regions):
            region_type = region.get('region_type', '')
            sequence_type = region.get('sequence_type', '')

            # Look for regions that are adapters or linkers
            if region_type in ['linker', 'primer', 'polyA', 'other'] or 'adapter' in region.get('name', '').lower():
                sequence = region.get('sequence', '')
                if sequence and sequence not in ['', 'N', 'X']:
                    adapters.append({
                        'name': region.get('name', f'adapter_{i}'),
                        'sequence': sequence,
                        'region_type': region_type,
                        'sequence_type': sequence_type,
                        'position': position,
                        'min_len': region.get('min_len', len(sequence)),
                        'max_len': region.get('max_len', len(sequence)),
                    })

            # Recursively search subregions
            if 'regions' in region and region['regions']:
                search_regions(region['regions'], position=f"{position}.{i}")

    # Search through all reads in the seqspec
    for read in seqspec.get('sequence_spec', []):
        if 'regions' in read:
            search_regions(read['regions'], position=read.get('read_id', 'unknown'))

    return adapters


def find_barcodes(seqspec: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Find all barcode/UMI regions in the seqspec.

    Returns:
        List of barcode information dictionaries
    """
    barcodes = []

    def search_regions(regions: List[Dict], read_id: str = "", position_in_read: int = 0):
        cumulative_pos = position_in_read
        for i, region in enumerate(regions):
            region_type = region.get('region_type', '')
            region_id = region.get('region_id', '')

            # Look for barcode or UMI regions
            if region_type in ['barcode', 'umi']:
                barcodes.append({
                    'name': region.get('name', region_id),
                    'region_id': region_id,
                    'region_type': region_type,
                    'read_id': read_id,
                    'position': cumulative_pos,
                    'min_len': region.get('min_len', 0),
                    'max_len': region.get('max_len', 0),
                    'onlist': region.get('onlist', None),
                    'sequence_type': region.get('sequence_type', ''),
                })

            # Update position for next region
            cumulative_pos += region.get('max_len', region.get('min_len', 0))

            # Recursively search subregions
            if 'regions' in region and region['regions']:
                search_regions(region['regions'], read_id, cumulative_pos)

    # Search through all reads in the seqspec
    for read in seqspec.get('sequence_spec', []):
        if 'regions' in read:
            search_regions(read['regions'], read_id=read.get('read_id', 'unknown'))

    return barcodes


def find_biological_sequence(seqspec: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Find the biological sequence region (typically cDNA or gDNA).

    Returns:
        Information about the biological sequence region
    """
    def search_regions(regions: List[Dict], read_id: str = "") -> Optional[Dict]:
        for region in regions:
            region_type = region.get('region_type', '')

            # Look for the actual biological sequence
            if region_type in ['cdna', 'gdna', 'rna', 'dna']:
                return {
                    'name': region.get('name', region_type),
                    'region_type': region_type,
                    'read_id': read_id,
                    'min_len': region.get('min_len', 0),
                    'max_len': region.get('max_len', 0),
                }

            # Recursively search subregions
            if 'regions' in region and region['regions']:
                result = search_regions(region['regions'], read_id)
                if result:
                    return result
        return None

    # Search through all reads in the seqspec
    for read in seqspec.get('sequence_spec', []):
        if 'regions' in read:
            result = search_regions(read['regions'], read_id=read.get('read_id', 'unknown'))
            if result:
                return result

    return None


def generate_adapter_config(seqspec: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate adapter configuration for cutadapt/fastp.

    Returns:
        Dictionary with adapter configuration
    """
    adapters = find_adapters(seqspec)

    config = {
        'adapters_5prime': [],
        'adapters_3prime': [],
        'adapters_anywhere': [],
    }

    for adapter in adapters:
        adapter_entry = {
            'name': adapter['name'],
            'sequence': adapter['sequence'],
        }

        # Try to determine adapter position based on region location
        position = adapter.get('position', '')
        if position.startswith('0') or 'first' in adapter['name'].lower() or '5' in adapter['name']:
            config['adapters_5prime'].append(adapter_entry)
        elif 'last' in position or 'last' in adapter['name'].lower() or '3' in adapter['name']:
            config['adapters_3prime'].append(adapter_entry)
        else:
            # If we can't determine, treat as "anywhere"
            config['adapters_anywhere'].append(adapter_entry)

    return config


def generate_processing_params(seqspec: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate processing parameters based on seqspec information.

    Returns:
        Dictionary with processing parameters
    """
    adapters = find_adapters(seqspec)
    barcodes = find_barcodes(seqspec)
    bio_seq = find_biological_sequence(seqspec)

    params = {
        'has_adapters': len(adapters) > 0,
        'adapter_count': len(adapters),
        'adapters': adapters,
        'has_barcodes': len(barcodes) > 0,
        'barcode_count': len(barcodes),
        'barcodes': barcodes,
        'biological_sequence': bio_seq,
        'requires_umi_extraction': any(bc['region_type'] == 'umi' for bc in barcodes),
        'requires_barcode_extraction': any(bc['region_type'] == 'barcode' for bc in barcodes),
    }

    # Determine minimum read length after trimming
    if bio_seq:
        params['min_read_length'] = bio_seq.get('min_len', 20)
    else:
        params['min_read_length'] = 20  # Default

    return params


def write_adapter_fasta(adapter_config: Dict[str, Any], output_path: str) -> int:
    """
    Write adapters to a FASTA file for fastp.

    Args:
        adapter_config: Dictionary with adapter lists
        output_path: Path to write the FASTA file

    Returns:
        Number of adapters written
    """
    adapter_count = 0

    with open(output_path, 'w') as f:
        # Write 5' adapters
        for adapter in adapter_config.get('adapters_5prime', []):
            seq = adapter['sequence']
            # Skip random sequences (UMIs, barcodes)
            if seq and not all(c in 'Nn' for c in seq):
                f.write(f">{adapter['name']}_5prime\n{seq}\n")
                adapter_count += 1

        # Write 3' adapters
        for adapter in adapter_config.get('adapters_3prime', []):
            seq = adapter['sequence']
            if seq and not all(c in 'Nn' for c in seq):
                f.write(f">{adapter['name']}_3prime\n{seq}\n")
                adapter_count += 1

        # Write "anywhere" adapters
        for adapter in adapter_config.get('adapters_anywhere', []):
            seq = adapter['sequence']
            if seq and not all(c in 'Nn' for c in seq):
                f.write(f">{adapter['name']}\n{seq}\n")
                adapter_count += 1

    return adapter_count


def main():
    parser = argparse.ArgumentParser(
        description='Process seqspec YAML files to extract adapter FASTA for fastp'
    )
    parser.add_argument(
        '--seqspec',
        required=True,
        help='Path to input seqspec YAML file'
    )
    parser.add_argument(
        '--output-fasta',
        required=True,
        help='Path to output adapter FASTA file'
    )

    args = parser.parse_args()

    try:
        # Parse seqspec
        seqspec = parse_seqspec(args.seqspec)

        # Generate adapter configuration
        adapter_config = generate_adapter_config(seqspec)

        # Write adapter FASTA
        adapter_count = write_adapter_fasta(adapter_config, args.output_fasta)

        print(f"Successfully processed seqspec file: {args.seqspec}")
        print(f"Wrote {adapter_count} adapters to {args.output_fasta}")
        print(f"  - {len(adapter_config['adapters_5prime'])} 5' adapters")
        print(f"  - {len(adapter_config['adapters_3prime'])} 3' adapters")
        print(f"  - {len(adapter_config['adapters_anywhere'])} other adapters")

    except Exception as e:
        print(f"Error processing seqspec: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
