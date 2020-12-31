from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# #individual sequences blast
from Bio import SeqIO

sequence_dict = SeqIO.to_dict(SeqIO.parse('fasta_sequences.fasta', 'fasta'))


# for name, sequence in sequence_dict.items():
#     blast_name = './Fasta/' + name +'.xml'
#     result_handle = NCBIWWW.qblast('blastn', 'nt', sequence.format("fasta"))
#     with open(blast_name, "w") as save_to:
#         save_to.write(result_handle.read())
            
        
#all sequences blast  
species = 'Treponema pallidum subsp. pallidum'  
query = 'nuccore pubmed[filter] NOT '+species+'[Organism]'    
sequence_data = open('fasta_sequences.fasta')
result_handle = NCBIWWW.qblast('blastn', 'nt', sequence_data.read(), entrez_query=query)
sequence_data.close()

with open("fasta_data.xml", "w") as save_to:
    save_to.write(result_handle.read())
    result_handle.close()

result_handle = open("fasta_data.xml", 'r')
blast_records = NCBIXML.parse(result_handle)
# blast_record = next(blast_records)
species = 'Treponema pallidum subsp. pallidum'

#Processing Blast Data 
blast_save = {}
for sequence in blast_records:
    E_VALUE_THRESH = 1
    unique = True
    for alignment in sequence.alignments:
        for hsp in alignment.hsps:
            title = str(alignment.title)
            if not hsp.expect < E_VALUE_THRESH:
                unique = False
    if unique:
       blast_save[sequence.query] = [alignment.title, hsp.query[0:23]]


# for title, sequence in blast_save.items():
#     if sequence == 'CCTGGTGTTGGTTACCGGCGTCG':  
#         print('yes')     
#         print('title')
##
# print(result_matches)
# print(fasta_seq, fasta_name)

# fasta_try = 'Example of a single sequence in FASTA/Pearson format:\n\n\n> sequence A\n' + fasta_try[1] + '\n\n'
# result_handle = NCBIWWW.qblast('blastn','nt', fasta_try)