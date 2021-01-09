from Bio.Blast import NCBIXML


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
    
    
    match_mismatch_number = {}
    for (query, matches) in blast_save.items():
        match_count = 0
        mismatch_count = 0
        for match in matches:
            (title, sequence) = match
            if species in title:
                match_count = match_count+1
            else:
                mismatch_count = mismatch_count + 1
        match_mismatch_number[query] = (match_count, mismatch_count)
     
    
    blast_results = {}
    for (query_match, matches) in blast_save.items():
        for (query_counts, counts) in match_mismatch_number.items():
            if query_match == query_counts:
                if query_match not in blast_results:
                    blast_results[query_match] = list()
                blast_results[query_match].append((matches, counts))
    return blast_results
                