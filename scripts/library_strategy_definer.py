import pandas as pd
import argparse

# Each dictionary stores the words for each library strategy and a score.
Ribo_seq = {'ribo-seq': 2, 'ribosome profil': 2, 'fp': 2, 'rpf': 2,  'ribolace': 2, 'riboseq': 2, 'ribosome protected fragments': 2, 'ribo': 2, 'footprint': 2, "harringtonine": 2,"lactimidomycin": 2, "initiating": 2,"pateamine" : 2,"harr": 2, "lac": 2, "cycloheximide": 2, "chx":2}
RNA_seq = {'rna-seq': 1, 'mrna': 1, 'rnaseq': 1, 'rna': 1}
Initiating_ribo = {"harringtonine": 1,"lactimidomycin": 1, "initiating": 1,"pateamine" : 1,"harr": 1, "lac": 1}
Elongating_ribo = {"Cyclohexamide": 1}

def get_scores_and_terms(field, terms_set):
    '''
    Given a field and a dictionary of terms, calculates the score for that specific field based on the terms found
    '''

    field_score = 0
    found_terms = []
    for term in terms_set:
        if term.casefold() in field.casefold():
            field_score += terms_set[term]
            if term not in found_terms:
                found_terms.append(term)
    
    return field_score, found_terms


def scores_evaluator(n, df):
    '''
    Using get_scores, determines from the title and eventually the description if the sample in ecam is Riboseq or RNA-seq.
    '''

    # Let's assign all the fields to a variable
    titl = df.at[n,"Title"].lower()
    des = df.at[n, "Description"].lower()

    #processing the info: I need to remove commas and dots -> technically not required if I use the "in" below. I look for substrings, not specific strings.

    ribo_titl_score, ribo_titl_found_terms = get_scores_and_terms(titl, Ribo_seq)
    RNA_titl_score, RNA_titl_found_terms = get_scores_and_terms(titl, RNA_seq)
    ribo_des_score = None
    RNA_des_score = None


    #  If there is a difference of less-than-one among the two scores, the description is also evaluated.
    if abs(ribo_titl_score - RNA_titl_score) < 1:

        ribo_des_score, ribo_des_found_terms = get_scores_and_terms(des, Ribo_seq)
        RNA_des_score, RNA_des_found_terms = get_scores_and_terms(des, RNA_seq)
        evidence = '[' + ','.join(ribo_titl_found_terms + ribo_des_found_terms + RNA_des_found_terms + RNA_titl_found_terms) + ']'

        if RNA_des_score - ribo_des_score > 0:
            return "RNA-seq study", evidence
        elif RNA_des_score - ribo_des_score < 0:
            return "Ribo-seq study", evidence
        elif RNA_des_score - ribo_des_score == 0:
            return "Inconclusive", evidence
    else:

        evidence = '[' + ','.join(ribo_titl_found_terms + RNA_titl_found_terms) + ']'
        if RNA_titl_score > ribo_titl_score:
            return "RNA-seq study", evidence
        elif RNA_titl_score < ribo_titl_score:
            return "Ribo-seq study", evidence
        else:
            return "Inconclusive", evidence


def define_ribosome_position(n,df):
    '''
    If the study is ribo-seq, checks for info about the ribosome position (stalling or elongating) using some keywords.
    '''
    prot = str(df.at[n,"Extraction_Protocol"]).lower()

    initiating_score = 0
    elongating_score = 0
    initiating_score, initiating_found_terms = get_scores_and_terms(prot, Initiating_ribo)
    elongating_score, elongating_found_terms = get_scores_and_terms(prot, Elongating_ribo)

    evidence = initiating_found_terms + elongating_found_terms

    if initiating_score > elongating_score:
        return "Ribo-seq study","Initiating", evidence
    elif initiating_score < elongating_score:
        return "Ribo-seq study","Elongating", evidence
    else:
        return "Ribo-seq study", "Inconclusive", "No Key Words Found"
    


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
        strategy_results = scores_evaluator(n, df)
        df.at[n,"Library_Strategy"] = strategy_results[0]
        df.at[n, "Library_Strategy_Evidence"] = strategy_results[1]

        if df.at[n,"Library_Strategy"] == "Ribo-seq study":
            ribosome_position_results = define_ribosome_position(n, df)
            df.at[n,"Library_Strategy"], df.at[n,"Ribosome_position"] = ribosome_position_results[0], ribosome_position_results[1]
            df.at[n, "Ribosome_Position_Evidence"] = ribosome_position_results[2]

        else: 
            df.at[n,"Ribosome_position"] = "NA"

    df.to_csv(file_path, index= False)
    



# However we must also identify which Ribo-Seq files are initiating and which are elongating ribosome protected fragments. 
# This is dictated by how the ribosomes are stalled. 

# Maybe this can be found in description or protocol?
ribosome_details_list = ["harringtonine","lactimidomycin", "initiating", "initiating ribosomes","pateamine","harr", "lac"]