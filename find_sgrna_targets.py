import data_blast
import search_pairs
import BLAST_pairs
import BLAST_pairs_taxid
import pandas as pd

date = "10012021"

#filenames & search species
chopchop_filename = "resultstrpI2.tsv"
fasta_filename = f"fasta_sequences_{date}.fasta"
fasta_filename_pairs = f"fasta_sequences_pairs{date}.fasta"
xml_filename = f"fasta_data_{date}.xml"
species, taxid=("Treponema pallidum subsp. pallidum", "161")

## search pairs
def search_different_pairs(chopchop_filename, fasta_filename, 
                           fasta_filename_pairs):
    header, data = search_pairs.read_chopchop_data(chopchop_filename)
    positive_grna_list, negative_grna_list = search_pairs.get_gRNA_data(header, data)
    result_matches, unique_matches_df = search_pairs.check_spacer_length_between_gRNAs(positive_grna_list, negative_grna_list)
    search_pairs.write_fasta_pairs_file(unique_matches_df, fasta_filename)
    search_pairs.write_fasta_pairs_file_for_pairs(result_matches, fasta_filename_pairs)
    return(result_matches)

# ## BLAST pairs
def blast_different_pairs(fasta_filename, taxid, xml_filename):
    BLAST_pairs_taxid.blast_all_sequences(fasta_filename, species, xml_filename)
    # blast_results = data_blast.process_blast_data(xml_filename, species)
    # return(blast_results)

# ## Merge pair & BLAST data
def get_pair_info(result_matches, blast_results):    
    total_results = pd.merge(result_matches,blast_results,left_on="match_no_positive", right_on="query")                                     
    total_results = total_results.drop(columns=['query'])                                             
    total_results = total_results.rename(columns={'match_count':'match_count_positive','mismatch_count':"mismatch_count_positive", "mismatch_titles":'mismatch_titles_positive', 'mismatch_sequences':'mismatch_sequences_positive'})
    total_results = pd.merge(total_results,blast_results,left_on="match_no_negative", right_on="query")                      
    total_results = total_results.rename(columns={'match_count':'match_count_negative','mismatch_count':"mismatch_count_negative", "mismatch_titles":'mismatch_titles_negative', 'mismatch_sequences':'mismatch_sequences_negative'})
    total_results.to_excel(f"{date}total_results.xlsx")
    
def master_code_that_runs_the_world(chopchop_filename, fasta_filename, xml_filename, species):
    result_matches = search_different_pairs(chopchop_filename, fasta_filename)
    blast_results = blast_different_pairs(fasta_filename, species, xml_filename)
    pairs_found = get_pair_info(result_matches, blast_results)
    return(pairs_found)
