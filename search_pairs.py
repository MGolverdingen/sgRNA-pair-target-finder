import pandas as pd
from Bio.Seq import Seq


def read_chopchop_data(chopchop_filename): 
    with open(chopchop_filename, 'r') as f:
        header = f.readline()
        data = f.readlines()
    return(header,data)
#%%

def get_gRNA_data(header:str, data:str, sequence_length:int):    
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
            position = position + sequence_length
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


#%%
def check_spacer_length_between_gRNAs(positive_g_rna_list:list,negative_g_rna_list:list, upper_limit:int, lower_limit:int):
    make_complement = lambda s: str(Seq(s).complement())
    UPPER_LIMIT = upper_limit
    LOWER_LIMIT = lower_limit
    # fasta_matches = ""
    result_matches =[]
    result_dict = {}
    count = 0
    for pos_g_rna in positive_g_rna_list:
        for neg_g_rna in negative_g_rna_list:
            distance = pos_g_rna['position'] - neg_g_rna['position']
            if distance <= UPPER_LIMIT and distance >= LOWER_LIMIT:
                # print(distance, pos_g_rna, neg_g_rna)
                # fasta_matches = fasta_matches + '>sequence A\n' + pos_g_rna['sequence'] + '\n\n>sequence B\n' + neg_g_rna['sequence']
                count = count + 1
                result_dict = {
                    "chopchop_quality_of_pair" : count,
                    "distance": distance,
                    "match_no_positive" : pos_g_rna['match_no'],
                    "match_no_negative" : neg_g_rna['match_no'],
                    "sequence_positive" : pos_g_rna['sequence'],
                    "sequence_negative" : neg_g_rna['sequence'],
                    "sequence_positive_c" : 
                        make_complement(pos_g_rna['sequence']), 
                    "sequence_negative_c" : 
                        make_complement(neg_g_rna['sequence'])
                }
                result_matches.append(result_dict)
                              
    matches_df = pd.DataFrame(result_matches)
    unique_matches_df = matches_df.drop_duplicates(['match_no_positive', 'match_no_negative'])        
    return (matches_df, unique_matches_df)

#%%
def write_fasta_pairs_file(unique_matches_df, fasta_filename:str):
  with open(fasta_filename, "w") as ofile:
      for index, row in unique_matches_df.iterrows():
          ofile.write(">" + row.match_no_positive +'\n' +
                      row.sequence_positive +'\n' + ">" +
                      row.match_no_negative + '_c' + '\n' +
                      row.sequence_negative_c + '\n')
                      # + ">" + row.match_no_positive +'_c' + '\n' +
                      # row.sequence_positive_c +'\n' + ">" +
                      # row.match_no_negative + '\n' +
                      # row.sequence_negative + '\n')

#%%   
def write_fasta_pairs_file_for_pairs(matches_df, fasta_filename_pairs:str):    
        with open(fasta_filename_pairs, "w") as ofile:
            for index, row in matches_df.iterrows():
                ofile.write(">" + row.match_no_positive + '-' +
                            row.match_no_negative + '_c' + '\n' +
                            row.sequence_positive +'-'+
                            row.sequence_negative_c + '\n')
                            # ">" + row.match_no_positive + '_c' + '-' +
                            # row.match_no_negative + '\n' +
                            # row.sequence_positive_c +'-'+
                            # row.sequence_negative + '\n')        
            
