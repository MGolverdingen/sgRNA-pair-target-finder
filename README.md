# sgRNA-pair-target-finder
## Description 
The sgRNA-pair-target finder is a custom python script that uses CHOPCHOP data of potential sgRNAs sequences to find pairs that have a PAM-out orientation and a spacer-length between a variable upper and lower limit.
By using NCBIWWW BLAST, the potential pairs are blasted to check their unqiueness for the species of interest. 
Using the CHOPCHOP RNA quality ranking and blast results, a ranking is made with the most suited sgRNAs for a paired CRISPR-system.
## How to Install
This tool works with Python v3. 
The required packages can be found in requirements.txt

## Example
chopchop_filename = "resultstrpI.tsv"  
fasta_filename = "fasta_sequences_pairs.fasta"  
xml_filename = "fasta_data.xml"  
species, taxid=("Treponema pallidum subsp. pallidum", "161")  
taxid_search = ["157", "9609"]  
sequence_length = 23  
upper_limit_spacer = 27  
lower_limit_spacer =17  

## Details
### Files
find_sgRNA_targets.py: 
file that runs other files and ensures an excel data output. 
*search_pairs.py*  
functions that read the .tsv file from CHOPCHOP, finds pairs with a given spacer length and a PAM-out orientation and gather the pairs in a .fasta file.     
*BLAST_pairs_taxid.py*  
functions that use NCBIWWW BLAST to search in specific genome databases for similar sequences as the sgRNA sequence. 
*datablast.py*  
function that uses BLAST data to count the number of mismatches of each pair.
### Input
chopchop_filename: a .tsv file with CHOPCHOP results 
fasta_filename: string with custom name for the .fasta file that is made in find_sgRNA_targets
xml_filename: string with custom name for the .xml file that is made in search_pairs.py 
species, taxid: tuple with strings of species and taxid following NCBI database nomenclature
taxid_search: list with strings of taxid of species you want to BLAST against
sequence_length: integer with sequence length of the single sgRNA
upper_limit_spacer: upperlimit of spacer length  
lower_limit_spacer: lowerlimit of spacer length
### Output
Excel Document