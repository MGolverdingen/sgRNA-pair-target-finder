import data_blast
import search_pairs
import BLAST_pairs

date = "09012021"

#filenames & search species
chopchop_filename = "resultstrpI2.tsv"
fasta_filename = f"fasta_sequences_{date}.fasta"
xml_filename = f"fasta_data_{date}.xml"
species="Treponema pallidum subsp. pallidum"



# ## search pairs
# header, data = search_pairs.read_chopchop_data(chopchop_filename)
# positive_grna_list, negative_grna_list = search_pairs.get_gRNA_data(header, data)
# result_matches, fasta_seq, fasta_name, unique_fasta_sequence, unique_fasta_name = search_pairs.check_spacer_length_between_gRNAs(positive_grna_list, negative_grna_list)
# search_pairs.write_fasta_pairs_file(unique_fasta_name, unique_fasta_sequence, fasta_filename)

## BLAST pairs

# BLAST_pairs.blast_all_sequences(fasta_filename, species, xml_filename)


blast_results = data_blast.process_blast_data(xml_filename, species)
print(blast_results)

