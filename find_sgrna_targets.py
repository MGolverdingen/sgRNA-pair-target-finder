import data_blast
import search_pairs
import BLAST_pairs_taxid
import pandas as pd
import numpy as np

date = "08022021"

#filenames & search species
chopchop_filename = "resultstrpI2.tsv"
fasta_filename = f"fasta_sequences_pairs{date}.fasta"
xml_filename = f"fasta_data_{date}.xml"
species, taxid=("Treponema pallidum subsp. pallidum", "161")
taxid_search = ["157", "9606"]
sequence_length = 23
upper_limit_spacer = 27
lower_limit_spacer =17

#%% search pairs
def search_different_pairs(chopchop_filename, 
                           fasta_filename_pairs, 
                           sequence_length,
                           upper_limit_spacer,
                           lower_limit_spacer):
    header, data = search_pairs.read_chopchop_data(chopchop_filename)
    positive_grna_list, negative_grna_list = search_pairs.get_gRNA_data(header, 
                                                                        data, 
                                                                        sequence_length)
    result_matches, unique_matches_df = search_pairs.check_spacer_length_between_gRNAs(positive_grna_list, 
                                                                                       negative_grna_list,
                                                                                       upper_limit_spacer,
                                                                                       lower_limit_spacer)
    search_pairs.write_fasta_pairs_file_for_pairs(result_matches, 
                                                  fasta_filename_pairs)

#%% BLAST pairs
def blast_different_pairs(fasta_filename, taxid, taxid_search, xml_filename, sequence_length):
    BLAST_pairs_taxid.blast_all_sequences(fasta_filename, taxid, taxid_search, xml_filename)
    blast_results = data_blast.process_blast_data(xml_filename, taxid, sequence_length)
    return(blast_results)

def find_full_length_seqs(in_list, full_seq_length):
    seq_lengths = list(map(len, in_list))
    full_length_seqs = (np.array(seq_lengths) == full_seq_length)
    return full_length_seqs
#%% Seperature full and partial matches for comparison of the pairs
def separate_full_and_partial_matches(df: pd.DataFrame, full_seq_length: int) -> pd.DataFrame:
    df['is_full_seq_length'] = df['mismatch_sequences'].map(lambda seq: find_full_length_seqs(seq, full_seq_length))
    df['full_length_mismatches'] = df.apply(lambda row: list(np.array(row.mismatch_sequences)[row.is_full_seq_length]), axis=1)
    df['full_length_titles'] = df.apply(lambda row: list(np.array(row.mismatch_titles)[row.is_full_seq_length]), axis=1)
    df['partial_length_mismatches'] = df.apply(lambda row: list(np.array(row.mismatch_sequences)[np.invert(row.is_full_seq_length)]), axis=1)
    df['partial_length_titles'] = df.apply(lambda row: list(np.array(row.mismatch_titles)[np.invert(row.is_full_seq_length)]), axis=1)
    df = df.drop("is_full_seq_length", axis=1)
    return df

#%% Write resuls to file
def write_results(blast_results):       
    blast_results.to_csv(f"{date}total_results.csv", index=False)
    
#%%
    
result_matches = search_different_pairs(chopchop_filename, fasta_filename, sequence_length, upper_limit_spacer, lower_limit_spacer)

#%%

blast_results = blast_different_pairs(fasta_filename, taxid, taxid_search, xml_filename, sequence_length)

#%%

separated_matches = separate_full_and_partial_matches(blast_results, sequence_length)

#%%

write_results(separated_matches)
