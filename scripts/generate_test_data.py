import random
from pathlib import Path
import click
from typing import List, Optional, Dict, Tuple
from dataclasses import dataclass
import json
from loguru import logger

@dataclass
class SequenceConfig:
    """Configuration for a sequence pattern (UMI or barcode)."""
    pattern_type: str  # 'UMI' or 'BARCODE'
    length: int
    position: str  # '5_prime' or '3_prime'
    sequence: Optional[str] = None  # For barcodes

class TestDataGenerator:
    def __init__(self, seed: Optional[int] = None):
        """Initialize generator with optional seed for reproducibility."""
        self.rng = random.Random(seed)
        
        # Define possible barcodes (can be expanded)
        self.possible_barcodes = {
            4: ['ACGT', 'TGCA', 'GATC', 'CTAG'],
            6: ['ACGTAC', 'TGCATG', 'GATCGA', 'CTAGCT'],
            8: ['ACGTACGT', 'TGCATGCA', 'GATCGATC', 'CTAGCTAG']
        }
    
    def generate_umi(self, length: int) -> str:
        """Generate a random UMI of specified length."""
        return ''.join(self.rng.choice('ACGT') for _ in range(length))
    
    def generate_read(self, length: int = 50) -> str:
        """Generate a random read of specified length."""
        return ''.join(self.rng.choice('ACGT') for _ in range(length))
    
    def add_sequence_pattern(self, read: str, pattern: SequenceConfig) -> str:
        """Add a UMI or barcode to a read at specified position."""
        if pattern.pattern_type == 'UMI':
            sequence = self.generate_umi(pattern.length)
        else:  # BARCODE
            sequence = self.rng.choice(self.possible_barcodes[pattern.length])
            
        if pattern.position == '5_prime':
            return sequence + read
        else:  # 3_prime
            return read + sequence
    
    def generate_quality_scores(self, length: int, min_qual: int = 30) -> str:
        """Generate quality scores for a read."""
        return ''.join(chr(self.rng.randint(min_qual, 40) + 33) for _ in range(length))
    
    def generate_fastq_record(self, read: str, read_num: int) -> List[str]:
        """Generate a complete FASTQ record."""
        return [
            f"@Read_{read_num}",
            read,
            "+",
            self.generate_quality_scores(len(read))
        ]

def generate_test_configs() -> List[Dict]:
    """Generate different test configurations."""
    configs = []
    
    # UMI configurations
    umi_lengths = [4, 6, 8]
    positions = ['5_prime', '3_prime']
    
    # Single UMI configurations
    for length in umi_lengths:
        for pos in positions:
            configs.append({
                'name': f'umi_{length}_{pos}',
                'patterns': [
                    SequenceConfig('UMI', length, pos)
                ]
            })
    
    # Dual UMI configurations
    for length1 in umi_lengths:
        for length2 in umi_lengths:
            configs.append({
                'name': f'dual_umi_{length1}_{length2}',
                'patterns': [
                    SequenceConfig('UMI', length1, '5_prime'),
                    SequenceConfig('UMI', length2, '3_prime')
                ]
            })
    
    # Barcode configurations
    barcode_lengths = [4, 6, 8]
    for length in barcode_lengths:
        for pos in positions:
            configs.append({
                'name': f'barcode_{length}_{pos}',
                'patterns': [
                    SequenceConfig('BARCODE', length, pos)
                ]
            })
    
    # Mixed UMI and barcode configurations
    for umi_len in umi_lengths:
        for barcode_len in barcode_lengths:
            configs.append({
                'name': f'umi_{umi_len}_5prime_barcode_{barcode_len}_3prime',
                'patterns': [
                    SequenceConfig('UMI', umi_len, '5_prime'),
                    SequenceConfig('BARCODE', barcode_len, '3_prime')
                ]
            })
    
    return configs

@click.command()
@click.option('--output-dir', '-o', required=True, type=click.Path(),
              help='Output directory for test FASTQ files')
@click.option('--num-reads', '-n', default=10000,
              help='Number of reads per test file')
@click.option('--read-length', '-l', default=30,
              help='Base read length (before adding UMIs/barcodes)')
@click.option('--seed', '-s', default=None, type=int,
              help='Random seed for reproducibility')
def main(outdir: str, num_reads: int, read_length: int, seed: Optional[int]):
    """Generate test FASTQ files with various UMI and barcode configurations."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Generating test files in {outdir}")
    
    # Initialize generator
    generator = TestDataGenerator(seed)
    
    # Get test configurations
    configs = generate_test_configs()
    
    # Generate metadata about test files
    metadata = {
        'parameters': {
            'num_reads': num_reads,
            'read_length': read_length,
            'seed': seed
        },
        'test_files': {}
    }
    
    # Generate test files for each configuration
    for config in configs:
        filename = f"test_{config['name']}.fastq"
        output_path = outdir / filename
        
        logger.info(f"Generating {filename}")
        
        with open(output_path, 'w') as f:
            for read_num in range(num_reads):
                # Generate base read
                read = generator.generate_read(read_length)
                
                # Add patterns according to configuration
                for pattern in config['patterns']:
                    read = generator.add_sequence_pattern(read, pattern)
                
                # Write FASTQ record
                fastq_record = generator.generate_fastq_record(read, read_num)
                f.write('\n'.join(fastq_record) + '\n')
        
        # Store configuration metadata
        metadata['test_files'][filename] = {
            'patterns': [
                {
                    'type': p.pattern_type,
                    'length': p.length,
                    'position': p.position
                }
                for p in config['patterns']
            ]
        }
    
    # Write metadata
    with open(outdir / 'test_metadata.json', 'w') as f:
        json.dump(metadata, f, indent=2)
    
    logger.info(f"Generated {len(configs)} test files")
    logger.info(f"Metadata written to {outdir}/test_metadata.json")

if __name__ == '__main__':
    main()