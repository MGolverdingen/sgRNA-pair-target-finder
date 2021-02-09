from Bio.Blast import NCBIWWW
from Bio import SeqIO

# #blast all sequences
# def individual_sequence_blast(filename, species):
#     sequence_dict = SeqIO.to_dict(SeqIO.parse(filename, 'fasta'))
       
#     for name, sequence in sequence_dict.items():
#         blast_name = './Fasta/' + name +'.xml'
#         result_handle = NCBIWWW.qblast('blastn', 'nt', sequence.format("fasta"))
#         with open(blast_name, "w") as save_to:
#             save_to.write(result_handle.read())
            
        
#blast specific taxids
def blast_all_sequences(fasta_filename:str, taxid:str, taxid_search:list, xml_filename:str):
    if taxid_search:
        taxid_query = '' 
        for i, taxid in enumerate(taxid_search):
            if i == 0:
                taxid_query = taxid_query + 'txid' + str(taxid) + '[ORGN]'
            else:
                taxid_query = taxid_query + ' OR ' + 'txid' + str(taxid) + '[ORGN]'
    sequence_data = open(fasta_filename)
    taxid_blast = ' all [filter] ('+ taxid_query+') NOT(environmental samples[ORGN] OR metagenomes[ORGN] OR txid'+str(taxid)+'[ORGN])  '
    result_handle = NCBIWWW.qblast('blastn', 'nt', sequence_data.read(), entrez_query=taxid_blast)
    sequence_data.close()
    
    with open(xml_filename, "w") as save_to:
        save_to.write(result_handle.read())
        result_handle.close()