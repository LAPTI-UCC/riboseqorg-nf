import pandas as pd
import argparse

# Each dictionary stores the words for each library strategy and a score.
Ribo_seq = {'ribo-seq': 1, 'ribosome profil': 1, 'fp': 1, 'rpf': 1,  'ribolace': 1, 'riboseq': 1, 'ribosome protected fragments': 1, 'ribo': 1, 'footprint': 1}
RNA_seq = {'rna-seq': 1, 'mrna': 1, 'rnaseq': 1, 'rna': 1}
Initiating_ribo = {"harringtonine": 1,"lactimidomycin": 1, "initiating": 1,"pateamine" : 1,"harr": 1, "lac": 1}
Elongating_ribo = {"Cyclohexamide": 1}

def get_score(field, terms_set):
    '''
    Given a field and a dictionary of terms, calculates the score for that specific field based on the terms found
    '''

    field_score = 0

    for term in terms_set:
        if term in field:
            #print("found term: ", term)
            field_score += terms_set[term]
            #print(field_score)
    
    return(field_score)


def scores_evaluator(n, df):
    '''
    Using get_scores, determines from the title and eventually the description if the sample in ecam is Riboseq or RNA-seq.
    '''

    # Let's assign all the fields to a variable
    titl = df.at[n,"Title"].lower()
    des = df.at[n, "Description"].lower()

    #processing the info: I need to remove commas and dots -> technically not required if I use the "in" below. I look for substrings, not specific strings.

    ribo_titl_score = get_score(titl, Ribo_seq)
    RNA_titl_score = get_score(titl, RNA_seq)

    ribo_des_score = None
    RNA_des_score = None

    #  If there is a difference of less-than-one among the two scores, the description is also evaluated.
    if abs(ribo_titl_score - RNA_titl_score) < 1:

        ribo_des_score = get_score(des, Ribo_seq)
        RNA_des_score = get_score(des, RNA_seq)

        if RNA_des_score - ribo_des_score > 0:
            return("RNA-seq study")
        elif RNA_des_score - ribo_des_score < 0:
            return("Ribo-seq study")
        elif RNA_des_score - ribo_des_score == 0:
            return("")

    if RNA_titl_score > ribo_titl_score:
        return("RNA-seq study")
    elif RNA_titl_score < ribo_titl_score:
        return("Ribo-seq study")

def define_ribosome_position(n,df):
    '''
    If the study is ribo-seq, checks for info about the ribosome position (stalling or elongating) using some keywords.
    '''
    prot = df.at[n,"Extraction_Protocol"].lower()

    initiating_score = 0
    elongating_score = 0
    initiating_score = get_score(prot, Initiating_ribo)
    elongating_score = get_score(prot, Elongating_ribo)

    if initiating_score > elongating_score:
        return ("Ribo-seq study","Initiating")
    elif initiating_score <= elongating_score:
        return ("Ribo-seq study","Elongating")


# function so that the new column is considered

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Retrieves the Library Strategy for a sepcific study, based on keywords")
    parser.add_argument(
        "path_to_csv_family",
        type=str,
        help="The path to the GSE_family.csv",
    )
    args = parser.parse_args()

    file_path = args.path_to_csv_family
    df = pd.read_csv(file_path)

    for n in df.index:
        df.at[n,"Library_Strategy"] = scores_evaluator(n, df)
        if df.at[n,"Library_Strategy"] == "Ribo-seq study":
            df.at[n,"Library_Strategy"], df.at[n,"Ribosome_position"] = define_ribosome_position(n, df)

    df.to_csv(file_path, index= False)
    



# However we must also identify which Ribo-Seq files are initiating and which are elongating ribosome protected fragments. 
# This is dictated by how the ribosomes are stalled. 

# Maybe this can be found in description or protocol?
ribosome_details_list = ["harringtonine","lactimidomycin", "initiating", "initiating ribosomes","pateamine","harr", "lac"]