from Bio.Blast import NCBIWWW
from Bio import SeqIO

def individual_sequence_blast(filename, species):
    sequence_dict = SeqIO.to_dict(SeqIO.parse(filename, 'fasta'))
       
    for name, sequence in sequence_dict.items():
        blast_name = './Fasta/' + name +'.xml'
        result_handle = NCBIWWW.qblast('blastn', 'nt', sequence.format("fasta"))
        with open(blast_name, "w") as save_to:
            save_to.write(result_handle.read())
            
        
#all sequences blast  
def blast_all_sequences(fasta_filename, taxid, xml_filename):
    # query = 'nuccore pubmed[filter] NOT '+species+'[Organism]'
    # query = '(none)'    
    sequence_data = open(fasta_filename)
    taxid_blast = ' all [filter] (txid157[ORGN] OR txid9606[ORGN]) NOT(environmental samples[ORGN] OR metagenomes[ORGN] OR txid161[ORGN])  '
    result_handle = NCBIWWW.qblast('blastn', 'nt', sequence_data.read(), entrez_query=taxid_blast)
    sequence_data.close()
    
    with open(xml_filename, "w") as save_to:
        save_to.write(result_handle.read())
        result_handle.close()