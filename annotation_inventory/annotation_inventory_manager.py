'''
This script facilitates all interactions with the annotation inventory. 
To date these processes include:
- Add organism 
- Remove Organism 
- Update Organism 

It is assumed the sqlite file is house in the same directoryas the script
and that it is named annotation_inventory.sqlite however, this can be specified. 
'''

import argparse
from operator import add
import sqlite3 
import os
# from calculating_chrom_sizes import get_chrom_sizes


def get_parser():
    '''
    return the argparse parser for this script 
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument("--operation", choices=['add', 'remove', 'update', 'set_primary_organism', 'add_gwips'])

    parser.add_argument("-o", help="Organism name. This is used as dictionary key for future reference", type=str) 
    parser.add_argument("-r", help="Path to rRNA bowtie indexes", type=str) 
    parser.add_argument("-t", help="Path to transcritome bowtie indexes", type=str) 
    parser.add_argument("-g", help="Path to genome bowtie indexes", type=str) 
    parser.add_argument("-f", help="Path to genome fasta file", type=str) 
    parser.add_argument("-s", help="Path to annotation inventory sqlite (default .annotation_inventory/annotation_inventory.sqlite)", type=str, default="annotation_inventory/annotation_inventory.sqlite") 
    parser.add_argument("-c", help="Path to chromosome sizes corresponding to genome fasta (-g)", type=str) 

    parser.add_argument("--scientific", help="Scientific name of the organism you want to set. Goes with '--operation' ", type=str)
    parser.add_argument("--gwips_db", help="Name of the gwips database for a given organism. Goes with '--operation' ", type=str)

    return parser


def get_db_cursor(connection):
    '''
    return a sqlite cursor for working with the db at the provided path
    '''
    cursor = connection.cursor()
    return cursor


def is_organism_in_use(organism, cursor):
    '''
    Return a boolean response describing whether the organism name is unique or not 
    '''
    organisms = cursor.execute("SELECT organism FROM annotation_inventory;").fetchall()
    return organism in [i[0] for i in organisms]


def handle_chrom_sizes(fasta_path):
    '''
    Create the missing chromosome sizes and return the path to that file 
    '''
    if not os.path.isabs(args.f):
        dir_path = os.path.dirname(os.path.abspath(args.f))
        chrom_sizes_path = dir_path + f"{args.f}.chrom_sizes"

    else:
        chrom_sizes_path = f"{args.f}.chrom_sizes"
                
    if not os.path.exists(chrom_sizes_path):
        get_chrom_sizes(args.f, chrom_sizes_path)
    return chrom_sizes_path


def add_organism(organism, rRNA_index, transcriptome_index, genome_index, genome_fasta, annotation_sqlite, chrom_sizes_file, db="annotation_inventory/annotation_inventory.sqlite"):
    '''
    Insert row with provided details into annotation inventory.
    '''
    connection = sqlite3.connect(db)
    cursor = get_db_cursor(connection)
    if not is_organism_in_use(organism, cursor):
        cursor.execute(
            f"INSERT INTO annotation_inventory VALUES('{organism}', '{rRNA_index}', '{transcriptome_index}', '{genome_index}', '{genome_fasta}', '{annotation_sqlite}', '{chrom_sizes_file}');"
            )
        connection.commit()
        connection.close()

    else:
        raise Exception("This organism already exists in the database. Use existing entry or create a new one with a more specific name")


def remove_organism(organism, db="annotation_inventory/annotation_inventory.sqlite"):
    '''
    Remove entry where organism = organism
    '''
    connection = sqlite3.connect(db)
    cursor = get_db_cursor(connection)
    if is_organism_in_use(organism, cursor):
        cursor.execute(
            f"DELETE FROM annotation_inventory WHERE organism = {organism}"
            )
        connection.commit()
        connection.close()

    else:
        results = cursor.execute("SELECT organism FROM annotation_inventory;").fetchall()
        raise Exception(f"{organism} isn't in the inventory. Maybe it is spelled differently?\n Choose from {results}")


def update_organism(organism, db="annotation_inventory/annotation_inventory.sqlite", **kwargs):
    '''
    Update the entry for the given organism using the parameters provided
    '''
    connection = sqlite3.connect(db)
    cursor = get_db_cursor(connection)
    if is_organism_in_use(organism, cursor):

        set_statement = "SET "
        for key, value in kwargs.items():
            if key == 'r':
                set_statement += f"rRNA_index='{value}',"
            elif key == 't':
                set_statement += f"transcriptome_index='{value}',"
            elif key == 'g':
                set_statement += f"genome_index='{value}',"
            elif key == 'f':
                set_statement += f"genome_fasta='{value}',"
            elif key == 's':
                set_statement += f"annotation_sqlite='{value}',"
            elif key == 'c':
                set_statement += f"chrom_sizes_file='{value}',"

        set_statement = ''.join(set_statement)
        
        cursor.execute(
            f"UPDATE annotation_inventory {set_statement} WHERE organism=='{organism}'"
            )
        connection.commit()
        connection.close()

    else:
        results = cursor.execute("SELECT organism FROM annotation_inventory;").fetchall()
        raise Exception(f"{organism} isn't in the inventory. Maybe it is spelled differently?\n Choose from {results}")


def set_primary_organism(organism, scientific, db='annotation_inventory/annotation_inventory.sqlite'):
    '''
    Set the primary organism for the given scientific name. (ie. 'Homo sapiens' primary organism 'human_hg19_gencode25')
    '''
    
    connection = sqlite3.connect(db)
    cursor = get_db_cursor(connection)

    organism_results = cursor.execute(f"SELECT COUNT(*) FROM annotation_inventory where organism='{organism}'").fetchall()[0][0]
    if organism_results < 1:
        raise Exception("This organism is not in the annotation inventory. Please add (using --operation add) first")

    scientific_results = cursor.execute(f"SELECT COUNT(*) FROM primary_organism where scientific_name='{scientific}'").fetchall()[0][0]
    if scientific_results < 1:
        cursor.execute(f"INSERT INTO primary_organism VALUES('{scientific}', '{organism}') ")
    else:
        cursor.execute(f"UPDATE primary_organism SET organism='{organism}")

    connection.commit()
    connection.close()


def add_gwips_organism(gwips_db, scientific, db='annotation_inventory/annotation_inventory.sqlite'):
    '''
    Adds the organism to the gwips_organism table

    Parameters:
        gwips_db: the name of the gwips database
        scientific: the scientific name of the organism
    '''
    connection = sqlite3.connect(db)
    cursor = get_db_cursor(connection)

    scientific_results = cursor.execute(f"SELECT COUNT(*) FROM gwips_organism where scientific_name='{scientific}'").fetchall()[0][0]

    if scientific_results < 1:
        cursor.execute(f"INSERT INTO gwips_organism VALUES('{scientific}', '{gwips_db}')")
    else:
        cursor.execute(f"UPDATE gwips_organism SET organism='{gwips_db}'")


def check_required_args(args, required: list) -> bool:
    '''
    Check that all required arguments are present
    
    Parameters:
        args: the arguments provided
        required: the required arguments
    '''
    for i in required:
        if not vars(args)[i]:
            raise Exception(f"Missing required argument: -{i}")
    return True

def run(args, db="annotation_inventory/annotation_inventory.sqlite"):
    '''
    Run the operation given the provided commands
    '''
    if args.operation == 'add':
        required = ['o', 'r','t','g','f']
        check_required_args(args, required)

        if not args.c:
            args.c = handle_chrom_sizes(args.f)

        add_organism(
            organism=args.o,
            rRNA_index=args.r,
            transcriptome_index=args.t,
            genome_index=args.g,
            genome_fasta=args.f,
            annotation_sqlite=args.s,
            chrom_sizes_file=args.c, 
            db=db
        )
    
    elif args.operation == "remove":
        required = ['o']
        check_required_args(args, required)

        remove_organism(organism=args.o, db=db)


    elif args.operation == "update":
        required = ['o']
        check_required_args(args, required)
        
        items_to_update = {key:value for key, value in vars(args).items() if value}

        for item in ['operation', 'o', 's']:
            items_to_update.pop(item)
            
        update_organism(organism=args.o, db=db, **items_to_update)
    
    elif args.operation == 'set_primary_organism':
        required = ['o', 'scientific']

        check_required_args(args, required)
        set_primary_organism(organism=args.o, db=db, scientific=args.scientific)

    elif args.operation == 'add_gwips':
        required = ['scientific', 'gwips_db']

        check_required_args(args, required)
        add_gwips_organism(gwips_db=args.gwips_db, scientific=args.scientific, db=db)
    
    else:
        raise Exception(f"Invalid operation: {args.operation}")


def run_prompted(args, db='annotation_inventory/annotation_inventory.sqlite'):
    '''
    Walk the user through parameter input and then run the program
    '''
    print("-"*90)
    args.db = input("What is the path to the annotation inventory database? (default: annotation_inventory/annotation_inventory.sqlite): ").strip(' ')
    print("-"*90)
    args.operation = input("Which operation would you like to perform on the database? ('add', 'remove' or 'update'): ").strip(' ')
    if args.operation not in ['add', 'remove', 'update', 'set_primary_organism', 'add_gwips']:
        raise Exception("invalid operation. Must be ('add', 'remove', 'update' or 'set_primary_organism')")
    
    print("-"*90)
    if args.operation == 'add':
        required = ['o', 'r','t','g','f', 's']
        options = ['c']

        args.o = input("What is the name of the organism you want to add? Be specific: ").strip(' ')
        print("-"*90)
        args.r = input("What is the path to the bowtie rRNA indices?: ").strip(' ')
        print("-"*90)
        args.t = input("What is the path to the bowtie transcriptome indices?: ").strip(' ')
        print("-"*90)
        args.g = input("What is the path to the bowtie genome indices?: ").strip(' ')
        print("-"*90)
        args.f = input("What is the path to the genome fasta file?: ").strip(' ')
        print("-"*90)
        args.s = input("What is the path to the trips annotation file?: ").strip(' ')



    elif args.operation == "remove":
        required = ['o']
        args.o = input("What is the name of the organism you want to add? Be specific: ").strip(' ')
        print("-"*90)


    elif args.operation == "update":
        required = ['o']
        options = ['r','t','g','f','c', 's']
        args.o = input("What is the name of the organism you want to update? Be specific: ").strip(' ')
        print("-"*90)
        args.r = input("What is the path to the bowtie rRNA indices?: ").strip(' ')
        print("-"*90)
        args.t = input("What is the path to the bowtie transcriptome indices?: ").strip(' ')
        print("-"*90)
        args.g = input("What is the path to the bowtie genome indices?: ").strip(' ')
        print("-"*90)
        args.f = input("What is the path to the genome fasta file?: ").strip(' ')
        print("-"*90)
        args.c = input("What is the path to the genomes chromosome sizes file?: ").strip(' ')
        print("-"*90)
        args.s = input("What is the path to the trips annotation file?: ").strip(' ')

        

    elif args.operation == "set_primary_organism":
        args.o = input("What is the name of the organism you want to set? Be specific: ").strip(' ')
        print("-"*90)
        args.scientific = input("What is the scientific name of the organism you want to set? eg. 'Homo sapiens': ").strip(' ')
        print("-"*90)

    run(args)


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if not args.operation:
        run_prompted(args)
    else:
        run(args)



'''
 python annotation_inventory/annotation_inventory_manager.py --operation update -o 'homo_sapiens' -r /mnt/data/indices/bowtie/ncRNA/homo_sapiens_hg38/homo_sapiens_hg38_rRNA -t /mnt/data/indices/bowtie/transcriptome/homo_sapiens_gencode39/homo_sapiens_gencode39_transcriptome -g /mnt/data/indices/bowtie/genome/homo_sapiens_hg38/homo_sapiens_hg38_genome -f /mnt/data/indices/bowtie/genome/homo_sapiens_hg38/hg38.fa -s /mnt/data/organism_sqlites/homo_sapiens_gencode39/homo_sapiens.gencode39.sqlite
'''