#!/usr/bin/env python3
"""
Extract offset information from ribosome profiling data stored in sqlitedict format.
"""

import argparse
import logging
from pathlib import Path
from typing import Dict
import pandas as pd
from sqlitedict import SqliteDict

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Extract offsets from ribosome profiling sqlitedict'
    )
    parser.add_argument(
        '-d', '--dict-path',
        required=True,
        help='Path to sqlitedict file'
    )
    parser.add_argument(
        '-o', '--outfile',
        type=Path,
        required=True,
        help='Output file path for results'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug logging'
    )
    parser.add_argument(
        '-e', '--end',
        default='fiveprime',
        choices=['fiveprime', 'threeprime'],
        help='End type to process (default: fiveprime)'
    )
    return parser.parse_args()


def main():
    """Main function"""
    args = get_args()
    
    # Set logging level
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    # Create parent directory if needed
    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        # Open sqlitedict and output file
        logger.info(f"Opening sqlitedict: {args.dict_path}")
        with SqliteDict(args.dict_path) as db, open(args.outfile, 'w') as out:
            # Write header
            out.write(f"read_length\toffset\n")
            
            # Extract and write offsets
            for read_length in db['offsets'][args.end]['offsets']:
                offset = db['offsets'][args.end]['offsets'][read_length]
                out.write(f"{read_length}\t{offset}\n")
                
        logger.info("Processing completed successfully")
            
    except Exception as e:
        logger.error(f"Error processing file: {e}")
        raise
if __name__ == "__main__":
    main()