import argparse
from Bio import Entrez

def get_bioproject_ids(gse):
    Entrez.email = "your_email@example.com"  # Replace with your email address
    handle = Entrez.esearch(db="gds", term=gse)
    record = Entrez.read(handle)
    handle.close()

    if record["IdList"]:
        gds_id = record["IdList"][0]
        handle = Entrez.esummary(db="gds", id=gds_id)
        summary_record = Entrez.read(handle)
        handle.close()
        srp = summary_record[0]["ExtRelations"][0]["TargetObject"]
        return srp
    else:
        print("No BioProject IDs found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve BioProject IDs associated with a GEO Accession (GSE)")
    parser.add_argument("gse", type=str, help="GEO Accession (GSE)")
    args = parser.parse_args()

    gse = args.gse
    bioproject_ids = get_bioproject_ids(gse)
    print(bioproject_ids)