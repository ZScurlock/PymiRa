#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PymiRa - Version 1
Created on Fri Mar 15 22:08:58 2024

@author: Zac Scurlock
"""
from collections import Counter
import numpy as np
from itertools import zip_longest, islice
from bisect import bisect_left, bisect_right
import gzip

def multi_open(file_path):
    return gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')

def upload_file(file_path):
    split_name = file_path.lower().split('.')
    if ('fasta' in split_name) or ('fa' in split_name):
        return parse_fasta(file_path)
    elif ('fastq' in split_name):
        return parse_fastq(file_path)
        
def parse_fasta(file_path):
    """
    Imports a FASTA / FASTA.gz file and substitutes U bases for T bases.
    """
    fasta_dict = {}
    sequence_lines = []
    read_name = None
    with multi_open(file_path) as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if sequence_lines:
                    fasta_dict[read_name] = "".join(sequence_lines).replace("U", "T")
                    sequence_lines = []
                read_name = line[1:]
            else:
                sequence_lines.append(line.replace("U", "T"))

        if sequence_lines:
            fasta_dict[read_name] = "".join(sequence_lines)

    return fasta_dict, list(fasta_dict.keys())



def parse_fastq(file_path):
    """
    Imports a FASTQ / FASTQ.gz file and substitutes U bases for T bases.
    """
    fastq_dict={}
    with multi_open(file_path) as file:
        while True:
            name = file.readline().rstrip()
            if not name:
                break
            seq = file.readline().rstrip().replace("U", "T")
            file.readline()
            file.readline()
            
            fastq_dict[name[1:]] = seq
    return fastq_dict, list(fastq_dict.keys())


def split_fasta_dict(fasta_dict, X):
    """
    Splits a FASTA dictionary into X smaller dictionaries.

    """
    if X < 1:
        raise ValueError("X must be at least 1")

    total_sequences = len(fasta_dict)
    split_size = total_sequences // X
    remainder = total_sequences % X

    splits = []

    read_names = list(fasta_dict.keys())
    start_index = 0

    for i in range(X):
        current_split_size = split_size + (1 if i < remainder else 0)

        split_dict = {
            read_names[j]: fasta_dict[read_names[j]]
            for j in range(start_index, start_index + current_split_size)
        }

        splits.append(split_dict)

        start_index += current_split_size

    return splits


def int_keys_best(l):
    seen = sorted(set(l))
    index = {v: i for i, v in enumerate(seen)}
    return [index[v] for v in l]

def bwt_from_suffix(string, s_array=None):
    if s_array is None:
        s_array = suffix_array2(string)
    return "".join(string[idx - 1] for idx in s_array)


def count_occurences(string, letters=None):
    count = 0
    result = {}

    c = Counter(string)
    if letters is None:
        letters = set(string)

    for letter in sorted(letters):
        result[letter] = count
        count += c[letter]
    return result


def update(begin, end, letter, lf_map, counts, string_length):
    beginning = counts[letter] + lf_map[letter][begin - 1] + 1
    ending = counts[letter] + lf_map[letter][end]
    return (beginning, ending)

def suffix_array2(s):
    n = len(s)
    k = 1
    line = np.array(int_keys_best(s), dtype='int64')    
    prev = np.zeros(n, dtype='int64')
    while line.max() < n - 1:
        prev[-k::] = -1
        prev[0:-k] = line[k::] + 1
        line = (n + 1) * line + prev
        k <<= 1
        _, line = np.unique(line, return_inverse=True)
    return line

def lf_mapping_2(bwt,letters=None):
    if letters is None:
        letters=set(bwt)
        
    n = len(bwt)
    results_array = np.zeros((n+2,len(letters)),dtype='int64')
    results_array_idx = {letter:i for i, letter in enumerate(letters)}
    
    counts_array=np.zeros(len(letters),dtype='int64')
    
    for idx,letter in enumerate(bwt):
        counts_array[results_array_idx[letter]] += 1
        results_array[idx] = counts_array
        
    results_array=results_array.T
    results = {letter:results_array[i] for i, letter in enumerate(letters)}
    
    return(results)

def generate_all(input_string, s_array=None, eos="$"):
    letters = set(input_string)
    counts = count_occurences(input_string, letters)

    input_string += eos
    
    if s_array is None:
        s_array = np.argsort(suffix_array2(input_string))
    bwt = bwt_from_suffix(input_string, s_array)
    lf_map = lf_mapping_2(bwt)

    return letters, bwt, lf_map, counts, s_array


def find_all(string, pattern):
    i = string.find(pattern)
    while i != -1:
        yield i
        i = string.find(pattern, i + 1)


def bwt_align(
    search_string,
    input_string,
    mismatches_5p=0,
    mismatches_3p=2,
    bwt_data=None,
    s_array=None,
):
    """
    Alignment function

    Parameters
    ----------
    search_string : str
        Read to align.
    input_string : str
        Reference sequence to align to
    mismatches_5p : TYPE, optional
        Number of mismatches permitted at the 5` end.
        (55% of read) The default is 0.
    mismatches_3p : TYPE, optional
        Number of mismatches permitted at the 3` end.
        (45% of read) The default is 2.
    bwt_data : list, optional
        BWT created data. The default is None.

    Returns
    -------
    list
        DESCRIPTION.

    """
    results = []
    result_new = []
    # inverted
    mismatch_flag = 1
    if len(search_string) == 0:
        return []
    if bwt_data is None:
        bwt_data = generate_all(input_string, s_array=s_array)

    letters, bwt, lf_map, count, s_array = bwt_data

    if len(letters) == 0:
        return []

    # If there are some letters not in letters, then disregard immediately
    if not set(search_string) <= letters:
        return []
    length = len(bwt)

    class Fuzzy(object):
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)

    fuz = [Fuzzy(search_string=search_string, begin=0, end=len(bwt) - 1, mismatches=0)]
    counter = 0
    while len(fuz) > 0:
        p = fuz.pop()
        counter = counter + 1
        searching = p.search_string[:-1]
        last = p.search_string[-1]

        # Possible letters the last one could be
        all_letters = [last] if p.mismatches == 0 else letters
        for letter in all_letters:

            begin, end = update(p.begin, p.end, letter, lf_map, count, length)

            if begin > end:
                mismatch_flag = 0
                if counter <= round(len(search_string) * 0.45):
                    search_string_75 = search_string[: round(len(search_string) * 0.55)]
                    search_string_25 = search_string[round(len(search_string) * 0.55) :]

                    # Searching 55% of read - 0 mismatches
                    fuz_5p = [
                        Fuzzy(
                            search_string=search_string_75,
                            begin=0,
                            end=len(bwt) - 1,
                            mismatches=mismatches_5p,
                        )
                    ]

                    while len(fuz_5p) > 0:
                        long_read = fuz_5p.pop()
                        search_long = long_read.search_string[:-1]

                        last_long = long_read.search_string[-1]

                        # Possible letters the last one could be
                        all_letters_long = (
                            [last_long] if long_read.mismatches == 0 else letters
                        )
                        for letter_long in all_letters_long:

                            begin_long, end_long = update(
                                long_read.begin,
                                long_read.end,
                                letter_long,
                                lf_map,
                                count,
                                length,
                            )
                            if begin_long > end_long:
                                return []

                            if begin_long <= end_long:
                                if len(search_long) == 0:
                                    results.extend(s_array[begin_long : end_long + 1])

                                else:
                                    miss = p.mismatches
                                    if letter_long != last_long:
                                        miss = max(0, long_read.mismatches - 1)
                                    fuz_5p.append(
                                        Fuzzy(
                                            search_string=search_long,
                                            begin=begin_long,
                                            end=end_long,
                                            mismatches=miss,
                                        )
                                    )

                    fuz_2 = [
                        Fuzzy(
                            search_string=search_string_25,
                            begin=0,
                            end=len(bwt) - 1,
                            mismatches=mismatches_3p,
                        )
                    ]
                    result_new = []
                    while len(fuz_2) > 0:
                        second = fuz_2.pop()
                        searching_2 = second.search_string[:-1]

                        last_2 = second.search_string[-1]
                        all_letters_2 = [last_2] if second.mismatches == 0 else letters

                        for base in all_letters_2:

                            begin_short, end_short = update(
                                second.begin, second.end, base, lf_map, count, length
                            )
                            if begin_short <= end_short:

                                if len(searching_2) == 0:

                                    result_new.extend(
                                        s_array[begin_short : end_short + 1]
                                    )

                                else:
                                    miss_2 = second.mismatches
                                    if base != last_2:
                                        miss_2 = max(0, second.mismatches - 1)
                                    fuz_2.append(
                                        Fuzzy(
                                            search_string=searching_2,
                                            begin=begin_short,
                                            end=end_short,
                                            mismatches=miss_2,
                                        )
                                    )

            if begin <= end:
                next
                if len(searching) != 0:
                    miss = p.mismatches
                    if letter != last:
                        miss = max(0, p.mismatches - 1)

                    fuz.append(
                        Fuzzy(
                            search_string=searching,
                            begin=begin,
                            end=end,
                            mismatches=miss,
                        )
                    )

                if len(searching) == 0:
                    results.extend(s_array[begin : end + 1])

    if results == []:
        results = None
    if result_new == []:
        result_new = None

    return results, result_new, mismatch_flag


class DecimalCounter(Counter):
    """
    A Class and method to calculate and assign counts to multiple
    mapping reads whereby 1 read maps to N loci and each loci receives
    1/N reads.
    """

    def update_division(self, iterable,num=1):
        increment = num / len(iterable)
        for elem in iterable:
            self[elem] += increment


def process_chunk(
    input_dict, ref_seq, ids_ref, bwt_data, integer_dict, mirna_flag=False, mismatches_5p=0, mismatches_3p=2, 
):
    """

    Parameters
    ----------
    input_dict : dict
        Dictionary of input_files, keys as read names and values as sequences
    ref_seq : str
        Reference sequence to align against
    ids_ref : list
        Reference sequence IDs
    bwt_data : list
        BWT created data
    integer_dict : dict
        Dictionary with sequences as keys and integers as values (number of times the sequence is found in the file)
    mirna_flag : bool
        Flag for to restrict miRNA alignments to either end of the hairpin sequence, adding either -5p / -3p notation. 
        The default is False.
    mismatches_5p : TYPE, optional
        DESCRIPTION. The default is 0.
    mismatches_3p : TYPE, optional
        DESCRIPTION. The default is 2.

    Returns
    -------
    test_dict : TYPE
        DESCRIPTION.
    res_dict : TYPE
        DESCRIPTION.

    """

    res_dict = {x: [] for x in input_dict}

    for i in input_dict:
        res_dict[i] = bwt_align(
            i,
            ref_seq,
            mismatches_5p=mismatches_5p,
            mismatches_3p=mismatches_3p,
            bwt_data=bwt_data,
        )


    # Removes reads with no hits
    res_dict = {k: v for k, v in res_dict.items() if any(v)}

    # Locates the position of all spaces in the reference
    arr = np.frombuffer(ref_seq.encode('ascii'), dtype=np.uint8)
    space_pos = np.where(arr == ord(' '))[0]
    space_pos +=1
    
    
    test_dict = DecimalCounter()
    # Refers the hit back to the reference sample
    for key in res_dict.copy():
        res = res_dict[key]
        if (len(res) > 0) & (len(res[0]) != 0):
            # Removes reads with 3 mismatches in main part of read.
            if (res[2] == 0) & (len(res[0]) >= 1) & (not bool(res[1])):
                #Unsuccess_dict.append(key) - Reason - too many mismatches in 5' end of read.
                res_dict.pop(key)
                continue
        else:
            res_dict.pop(key)
            #Unsuccess_dict.append(key) - Reason - should not be here
            continue

        # Mismatch filter is 0 as mismatch and 1 as no mismatch
        main_pos = np.array(res[0])
        pos_3p, mm = res[1:]

        # True if there was a mismatch and there is a match in the 3`.
        if bool(pos_3p):
            np_main_pos = np.sort(np.array(main_pos))
            np_pos_3p = np.sort(np.array(pos_3p))

            # Removes any 3` match less than the lowest 5` match - these can't be real.
            np_pos_3p = np_pos_3p[np_pos_3p >= np_main_pos.min()]
            if len(np_pos_3p) == 0:
                res_dict.pop(key)
                #Unsuccess_dict.append(key) - Reason - no valid 3p end alignment
                continue

            # Finds the value immediately greater than each main_pos value
            indices = [
                bisect_right(np_pos_3p, np_main_pos[x]) - 1
                if bisect_right(np_pos_3p, np_main_pos[x]) >= len(np_pos_3p)
                else bisect_right(np_pos_3p, np_main_pos[x])
                for x in range(len(np_main_pos))
            ]
            main_pos_2 = []
            # For each main_pos value, it takes the closest match at the 3` and looks if there is a space - if so - then it hasn't matched to the same reference seq
            for ind, val in enumerate(indices):
                if (ref_seq[np_main_pos[ind] : np_pos_3p[val]].find(" ") == -1) & (
                    len(ref_seq[np_main_pos[ind] : np_pos_3p[val]]) > 0
                ):
                    main_pos_2.append(np_main_pos[ind])
                else:
                    continue

            main_pos = np.array(main_pos_2)
            if len(main_pos) < 1:
                res_dict.pop(key)
                continue
        
        # 
        mask=np.isin(main_pos,space_pos)
        main_pos[mask] +=1
        
        read_ind = np.unique(np.searchsorted(space_pos,main_pos,side='left'))
        ids_ind = read_ind-1
        
        ids = [ids_ref[x] for x in ids_ind]
        subject_seq=[]
        
        for actual ,num in enumerate(read_ind):
            subject_seq.append(ref_seq[space_pos[num - 1] : space_pos[num] - 1])
            #if ref_seq[space_pos[num - 1] : space_pos[num] - 1] != ref_dict[ids[actual]]:
            #    print('NOT MATCHING', key, ref_ids[actual])
                
        
        #names = []
        int_name = [x for x in ids]
        if mirna_flag:
            orient=[]
            newlist=[]
            for i in range(len(ids)):
            
                #id_counter+=1
            #Remove matches in middle of miRNA hairpins
                min_bound = (len(subject_seq[i])/2 + space_pos[read_ind[i]-1]) - 4.5
                max_bound = (len(subject_seq[i])/2 + space_pos[read_ind[i]-1]) + 4.5
                if not (main_pos[i] > min_bound) or not (main_pos[i] < max_bound):
                    newlist.append(ids[i])
                    
                    #Add orientations
                    if main_pos[i] + 9 < space_pos[read_ind[i]-1] + (
                        len(subject_seq[i]) / 2
                    ):
                        orient.append("-5p")
                    else:
                        orient.append("-3p")
                        
            if len(newlist) == 0:
                res_dict.pop(key)
                continue
            
            int_name = [ids[n].split(" MI")[0] + orient[n] for n in range(len(newlist))]
        
        res_dict[key] = int_name
        #res_dict_counter+=1
        if len(int_name) > 1:
            test_dict.update_division(res_dict[key],integer_dict[key])
        else:
            test_dict.update(res_dict[key]*integer_dict[key])
        continue
    return test_dict,res_dict
