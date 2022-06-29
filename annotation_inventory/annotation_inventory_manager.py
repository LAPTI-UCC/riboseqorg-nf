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
from calculating_chrom_sizes import get_chrom_sizes


def get_parser():
    '''
    return the argparse parser for this script 
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument("--operation", choices=['add', 'remove', 'update', 'set_primary_organism'])

    parser.add_argument("-o", help="Organism name. This is used as dictionary key for future reference", type=str) 
    parser.add_argument("-r", help="Path to rRNA bowtie indexes", type=str) 
    parser.add_argument("-t", help="Path to transcritome bowtie indexes", type=str) 
    parser.add_argument("-g", help="Path to genome bowtie indexes", type=str) 
    parser.add_argument("-f", help="Path to genome fasta file", type=str) 
    parser.add_argument("-s", help="Path to annotation inventory sqlite (default .annotation_inventory.sqlite)", type=str, default="annotation_inventory.sqlite") 
    parser.add_argument("-c", help="Path to chromosome sizes corresponding to genome fasta (-g)", type=str) 

    parser.add_argument("--scientific", help="Scientific name of the organism you want to set. Goes with '--operation ", type=str)
  
    return parser


def get_db_curosor(connection):
    '''
    return a sqlite cursor for working with the db at the provided path
    '''
    cursor = connection.cursor()
    return cursor


def is_organism_in_use(organism, curosr):
    '''
    Return a boolean response describing whether the organism name is unique or not 
    '''
    organisms = curosr.execute("SELECT organism FROM annotation_inventory;").fetchall()
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


def add_organism(organism, rRNA_index, transcriptome_index, genome_index, genome_fasta, annotation_sqlite, chrom_sizes_file, db="annotation_inventory.sqlite"):
    '''
    Insert row with provided details into annotation inventory.
    '''
    connection = sqlite3.connect(db)
    cursor = get_db_curosor(connection)
    if not is_organism_in_use(organism, cursor):
        cursor.execute(
            f"INSERT INTO annotation_inventory VALUES('{organism}', '{rRNA_index}', '{transcriptome_index}', '{genome_index}', '{genome_fasta}', '{annotation_sqlite}', '{chrom_sizes_file}');"
            )
        connection.commit()
        connection.close()

    else:
        raise Exception("This organism already exists in the database. Use existing entry or create a new one with a more specific name")


def remove_organism(organism, db="annotation_inventory.sqlite"):
    '''
    Remove entry where organism = organism
    '''
    connection = sqlite3.connect(db)
    cursor = get_db_curosor(connection)
    if is_organism_in_use(organism, cursor):
        cursor.execute(
            f"DELETE FROM annotation_inventory WHERE organism = {organism}"
            )
        connection.commit()
        connection.close()

    else:
        raise Exception(f"{organism} isn't in the inventory. Maybe it is spelled differently?")


def update_organim(organism, db="annotation_inventory.sqlite", **kwargs):
    '''
    Update the entry for the given organims using the parameters provided
    '''
    connection = sqlite3.connect(db)
    cursor = get_db_curosor(connection)
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

        set_statement = ''.join(set_statement.split(','))
        
        cursor.execute(
            f"UPDATE annotation_inventory {set_statement} WHERE organism=='{organism}'"
            )
        connection.commit()
        connection.close()

    else:
        raise Exception(f"{organism} isn't in the inventory. Maybe it is spelled differently?")


def set_primary_organism(organism, scientific, db='annotation_inventory.sqlite'):
    '''
    Set the primary organism for the given scientific name. (ie. 'Homo sapiens' primary organism 'human_hg19_gencode25')
    '''
    
    connection = sqlite3.connect(db)
    cursor = get_db_curosor(connection)

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

def run(args):
    '''
    Run the operation given the provided commands
    '''
    if args.operation == 'add':
        required = ['o', 'r','t','g','f']
        for i in required:
            if not vars(args)[i]:
                raise Exception(f"Missing required argument: -{i}")

        if not args.c:
            args.c = handle_chrom_sizes(args.f)

        add_organism(
            organism=args.o,
            rRNA_index=args.r,
            transcriptome_index=args.t,
            genome_index=args.g,
            genome_fasta=args.f,
            annotation_sqlite=args.s,
            chrom_sizes_file=args.c
        )
    
    elif args.operation == "remove":
        required = ['o']
        for i in required:
            if not vars(args)[i]:
                raise Exception(f"Missing required argument: -{i}")

        remove_organism(organism=args.o, db=args.s)


    elif args.operation == "update":
        required = ['o']
        for i in required:
            if not vars(args)[i]:
                raise Exception(f"Missing required argument: -{i}")
        
        items_to_update = {key:value for key, value in vars(args).items() if value}

        for item in ['operation', 'o', 's']:
            items_to_update.pop(item)
            
        update_organim(organism=args.o, db=args.s, **items_to_update)
    
    elif args.operation == 'set_primary_organism':
        required = ['o', 'scientific']

        for i in required:
            if not vars(args)[i]:
                raise Exception(f"Missing required argument: -{i}")

        set_primary_organism(organism=args.o, db=args.s, scientific=args.scientific)
    
    else:
        raise Exception(f"Invalid operation: {args.operation}")


def run_prompted(args, db='annotation_inventory.sqlite'):
    '''
    Walk the user through parameter input and then run the program
    '''
    print("-"*90)
    args.operation = input("Which operation would you like to perform on the database? ('add', 'remove' or 'update'): ").strip(' ')
    if args.operation not in ['add', 'remove', 'update', 'set_primary_organim']:
        raise Exception("invalid opetion. Must be ('add', 'remove', 'update' or 'set_primary_organim')")
    
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

        

    elif args.operation == "set_primary_organim":
        args.o = input("What is the name of the organism you want to set? Be specific: ").strip(' ')
        print("-"*90)
        args.scientific = input("What is the scientific name of the organism you want to set? eg. 'Homo sapiens': ").strip(' ')
        print("-"*90)

    run(args)


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if not args.operation:
        run_prompted(args, db=args.s)
    else:
        run(args)



