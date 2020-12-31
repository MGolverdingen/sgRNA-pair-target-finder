with open('resultstrpI2.tsv', 'r') as f:
    header = f.readline()
    data = f.readlines()
##
def get_gRNA_data(header, data):    
    MATCH_NUMBER_FIELD = 0
    POSITION_FIELD = 2
    SIGN_FIELD = 3
    SEQUENCE_FIELD = 1
    g_rna_list = []
    
    for line in data:
        fields = line.split('\t')
        position_field = fields[POSITION_FIELD]
        match_number = fields[MATCH_NUMBER_FIELD]
        sign = fields[SIGN_FIELD]
        sequence = fields[SEQUENCE_FIELD]
        position = int(position_field.split(":")[1])
        if sign == '-':
            position = position + 23
        g_rna = {
            "match_no": match_number,
            "sign": sign,
            "position": position,
            "sequence": sequence
        }
        g_rna_list.append(g_rna)
    ##
    positive_g_rna_list = [g for g in g_rna_list if g['sign'] == '+']
    negative_g_rna_list = [g for g in g_rna_list if g['sign'] == '-']
    return(positive_g_rna_list, negative_g_rna_list)
##
def check_spacer_length_between_gRNAs(positive_g_rna_list,negative_g_rna_list):
    UPPER_LIMIT = 35
    LOWER_LIMIT = 5
    result_matches = []
    # fasta_matches = ""
    fasta_seq = []
    fasta_name = []
    for pos_g_rna in positive_g_rna_list:
        for neg_g_rna in negative_g_rna_list:
            distance = pos_g_rna['position'] - neg_g_rna['position']
            if distance <= UPPER_LIMIT and distance >= LOWER_LIMIT:
                # print(distance, pos_g_rna, neg_g_rna)
                # fasta_matches = fasta_matches + '>sequence A\n' + pos_g_rna['sequence'] + '\n\n>sequence B\n' + neg_g_rna['sequence']
                result_matches.append([distance, pos_g_rna,neg_g_rna])
                fasta_seq.append(pos_g_rna['sequence'])
                fasta_seq.append(neg_g_rna['sequence'])
                fasta_name.append(pos_g_rna['match_no'])
                fasta_name.append(neg_g_rna['match_no'])
  
    fasta_name_uqe = list(set(fasta_name))
    fasta_seq_uqe = list(set(fasta_seq))
    return (fasta_seq, fasta_name, fasta_name_uqe, fasta_seq_uqe)

with open("fasta_pairs.fasta", "w") as ofile:
    for i in range(0,len(fasta_seq),2):
        ofile.write(">" + fasta_name[i] + ' ' + fasta_name[i+1] + '\n' + fasta_seq[i] + '\n' + fasta_seq[i+1] + '\n')

for i in range(len(fasta_seq)):
    fasta_dir= 'D:/Google Drive/01 Studie/02 MEP/01 Literature/02 Syphilis/Python/Fasta/' + fasta_name[i]+ '.fasta'
    print(fasta_dir)
    with open(fasta_dir, "w") as file:
        file.write(">" + fasta_name[i] + '\n' + fasta_seq[i] + '\n')

with open("fasta_sequences.fasta", "w") as ofile:
    for i in range(len(fasta_seq_uqe)):
        ofile.write(">" + fasta_name_uqe[i] + '\n' + fasta_seq_uqe[i] + '\n')

##############################################
# sequence_data = open('fasta_pairs.fasta')
# result_handle = NCBIWWW.qblast('blastn', 'nt', sequence_data.read())
# sequence_data.close()


# with open("my_blast.xml", "w") as save_to:
#     save_to.write(result_handle.read())
#     result_handle.close()

# result_handle = open("my_blast.xml", 'r')
# blast_records = NCBIXML.parse(result_handle)
# blast_record = next(blast_records)


# blast_save = {}
# E_VALUE_THRESH = 1
# for alignment in blast_record.alignments:
#     for hsp in alignment.hsps:
#         if hsp.expect < E_VALUE_THRESH:
#             print('****Alignment****')
#             print('sequence:', alignment.title)
#             print('length:', alignment.length)
#             print('e value:', hsp.expect)
#             print(hsp.query[0:23] + '...')
#             print(hsp.match[0:23] + '...')
#             print(hsp.sbjct[0:23] + '...')
#             blast_save[alignment.title] = hsp.query[0:23]


# for title, sequence in blast_save.items():
#     if sequence == 'CCTGGTGTTGGTTACCGGCGTCG':  
#         print('yes')     
#         print('title')
##
# print(result_matches)
# print(fasta_seq, fasta_name)

# fasta_try = 'Example of a single sequence in FASTA/Pearson format:\n\n\n> sequence A\n' + fasta_try[1] + '\n\n'
# result_handle = NCBIWWW.qblast('blastn','nt', fasta_try)

