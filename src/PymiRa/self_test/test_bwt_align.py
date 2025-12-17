#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: zac
"""
from PymiRa.functions import generate_all, split_fasta_dict, process_chunk,merge_results
from collections import defaultdict

def test_end_to_end():
    ref = ' GGGACGTACGTACGTACGTACGTCCCGGGACGTACGTACGTACGTACGTTTT '
    ref_ids = ['test_ref']
    input_dict = {'test_exact': 'ACGTACGTACGTACGT',
             'test_1mm_3p': 'ACGTACGTACGTACGA',
             'test_2mm_3p': 'ACGTACGTACGTTTGT',
             'test_3mm_3p': 'ACGTACGTACGATCTT',
             'test_1mm_5p': 'ACGTGCGTACGTACGT',
             'test_1mm_5p_2mm_3p': 'ACTTACGTACGTTGGT'}
    
    seq_dict = defaultdict(list)
    for name,seq in input_dict.items():
        seq_dict[seq].append(name)
    
    int_dict = {k:len(v) for k,v in seq_dict.items()}
    test_cont = split_fasta_dict(int_dict,2)
    letters, bwt, lf_map, count, s_array = generate_all(ref)
    required = [letters, bwt, lf_map, count, s_array]
    
    cont = []
    cont.append(process_chunk(input_dict = test_cont[0],ref_seq = ref, ids_ref = ref_ids,
                                 bwt_data = required,integer_dict = int_dict, mirna_flag = False,
                                 mismatches_3p=2))
    cont.append(process_chunk(input_dict = test_cont[1],ref_seq = ref, ids_ref = ref_ids,
                                 bwt_data = required,integer_dict = int_dict, mirna_flag = False,
                                 mismatches_3p=2))
    
    final = merge_results(cont)
    log = {}
    for up in range(len(cont)):
        log.update(cont[up][1])
    
    new_log = {}
    for k,v in log.items():
        all_read = seq_dict[k]
        for read in all_read:
            new_log[read] = v
    
    if len(new_log.keys()) != 3 or sum(final.values()) != 3:
        return False, f"Unexpected end-to-end results:{len(new_log.keys())} should equal 3."
    
    return True, "Alignment test passed."

test_end_to_end()
