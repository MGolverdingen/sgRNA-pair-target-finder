from Bio.Blast import NCBIXML
import pandas as pd

def process_blast_data(xml_filename:str, species:str):
    result_handle = open(xml_filename, 'r')
    blast_records = NCBIXML.parse(result_handle)
    
    
    #Processing Blast Data 
    blast_save = {}
    for sequence in blast_records:
        E_VALUE_THRESH = 1
        unique = True
        for alignment in sequence.alignments:
           # print(str(sequence.query))
            for hsp in alignment.hsps:
                title = str(alignment.title)
                if not hsp.expect < E_VALUE_THRESH:
                    unique = False
                if unique:
                    if sequence.query not in blast_save:
                        blast_save[sequence.query] = list()
                    blast_save[sequence.query].append((alignment.title, hsp.query[0:23]))
    
    # need to add the number of matches in the blast, not sure if this BLAST results only 100% matches or also others
    
    match_mismatch_number = {}
    match_mismatch_list = []
    for (query, matches) in blast_save.items():
        match_count = 0
        mismatch_count = 0
        mismatch_titles = []
        mismatch_sequence = []
        for match in matches:
            (title, sequence) = match
            if species in title:
                match_count = match_count+1
            else:
                mismatch_count = mismatch_count + 1
                mismatch_titles.append(title)
                mismatch_sequence.append(sequence)
        match_mismatch_number = {"query" : query, 
                                 "match_count": match_count,
                                 "mismatch_count": mismatch_count,
                                 "mismatch_titles": mismatch_titles,
                                 "mismatch_sequences": mismatch_sequence
        }
        match_mismatch_list.append(match_mismatch_number)
    matches_df = pd.DataFrame(match_mismatch_list)
    return(matches_df)