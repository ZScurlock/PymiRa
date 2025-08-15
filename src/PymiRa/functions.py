#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PymiRa - Version 1
Created on Fri Mar 15 22:08:58 2024

@author: Zac Scurlock
"""
import multiprocessing
from collections import Counter
import numpy as np
import pandas as pd
from itertools import zip_longest, islice
import re
from bisect import bisect_left, bisect_right
import time
from functools import partial
import sys

def parse_fasta(file_path):
    fasta_dict={}
    sequence_lines = []
    read_name=None
    print(f'Importing Fasta file {file_path}')
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('>'):
                if sequence_lines:
                    fasta_dict[read_name] = ''.join(sequence_lines).replace('U', 'T')
                    sequence_lines = []
                read_name = line[1:] 
            else:
                sequence_lines.append(line)
                
        if sequence_lines:
            fasta_dict[read_name] = ''.join(sequence_lines)

    return fasta_dict, list(fasta_dict.keys())

def split_fasta_dict(fasta_dict, X):
    """
    Splits a FASTA dictionary into X smaller dictionaries.

    :param fasta_dict: The original FASTA dictionary with read names and sequences.
    :param X: The number of splits to create.
    :return: A list of dictionaries split into X parts.
    """
    if X < 1:
        raise ValueError("X must be at least 1")

    total_sequences = len(fasta_dict)
    split_size = total_sequences // X
    remainder = total_sequences % X

    splits = []

    read_names = list(fasta_dict.keys())  # Get the list of read names
    start_index = 0

    for i in range(X):
        # Determine the size for the current split
        current_split_size = split_size + (1 if i < remainder else 0)

        split_dict = {read_names[j]: fasta_dict[read_names[j]]
                      for j in range(start_index, start_index + current_split_size)}

        splits.append(split_dict)

        # Update the starting index 
        start_index += current_split_size

    return splits


def to_int_keys_best(l):

    seen = set()
    ls = []
    for e in l:
        if not e in seen:
            ls.append(e)
            seen.add(e)
    ls.sort()
    index = {v: i for i, v in enumerate(ls)}
    return [index[v] for v in l]

def to_int_keys_best2(l):
    seen=set(l)
    ls=list(seen)
    ls.sort()
    index={v: i for i, v in enumerate(ls)}
    return [index[v] for v in l]

def to_int_keys_best3(l):
    seen=sorted(set(l))
    index={v:i for i,v in enumerate(seen)}
    return [index[v] for v in l]

def suffix_array(s):

    n = len(s)
    k = 1
    line = to_int_keys_best3(s)
    while max(line) < n - 1:
        line = to_int_keys_best3(
            [a * (n + 1) + b + 1
             for (a, b) in
             zip_longest(line, islice(line, k, None),
                         fillvalue=-1)])
        k <<= 1
    return line

def suffix_array3(s):
    n=len(s)
    k=1
    line = np.array(to_int_keys_best3(s), dtype='int64')
    #prev = np.zeros(n)
    while line.max() <n-1:
        line = np.array(to_int_keys_best3(
            [a * (n + 1) + b + 1
             for (a, b) in
             zip_longest(line, islice(line, k, None),
                         fillvalue=-1)]))
        k<<=1
    return line

def suffix_array2(s):
    '''
    Memory needed for this solution is 2*2n*64/8/1000000 mb
    creates a super fast suffix array for a string in the following steps
    1) convert the string to a list of numbers with the 'to_int_keys_best_new' function
    2) since 'to_int_keys_best_new' returns the numbers [0, N] with N being the number of unique elements in the list,
       we know that if the max element in that list is the size-1, each element in the list is unique and sorted.
    3) iteratively update each element in the list of numbers using "a * (n + 1) + b + 1" where n is the size of 
       the string to be sorted, a is the current number in the list at index i and b is the number in the list at index i+k
       k is incremented according to 2**iteration. The idea behind this formula is that you try to find how many of the numbers
       in the list you need to compare, in order to be able to sort it. We can increment k with 2**iteration since the number at 
       each index becomes the stacked 'sorted value' which doubles every iteration.
       The easiest way to see this is as a 'lego tower' start with a sequence of blocks.
       See figure at the bottom of the page.
    4) the np.unique function does the same as 'to_int_keys_best', but it only works for numpy arrays, but is much(!!!) faster.
    5) if we want to improve for memory we can replace the prev and line bit with:
        line[0:-k] = (n + 1) * line[0:-k] + line[k::] + 1
        line[-k::] *= (n + 1)

        however, that's arguibly less readable
    '''
    n = len(s)
    k = 1
    line = np.array(to_int_keys_best2(s), dtype='int64')    
    prev = np.zeros(n, dtype='int64')
    while line.max() < n - 1:
        prev[-k::] = -1
        prev[0:-k] = line[k::] + 1
        line = (n + 1) * line + prev
        k <<= 1
        _, line = np.unique(line, return_inverse=True)
    return line

def inverse_array(l):
    n = len(l)
    ans = [0] * n
    for i in range(n):
        ans[l[i]] = i
    return ans

#np.zeros(len(input_string))

def bwt_from_suffix(string, s_array=None):
    if s_array is None:
        s_array = suffix_array3(string)
    return("".join(string[idx - 1] for idx in s_array))


def lf_mapping(bwt, letters=None):
    if letters is None:
        letters = set(bwt)
        
    result = {letter:[0] for letter in letters}
    result[bwt[0]] = [1]
    for letter in bwt[1:]:
        for i, j in result.items():
            j.append(j[-1] + (i == letter))
    return(result)



def count_occurences(string,letters=None):
    count = 0
    result = {}
        
    c = Counter(string)
    if letters is None:
        letters=set(string)
        
    for letter in sorted(letters):
        result[letter] = count
        count += c[letter]
    return result


def update(begin, end, letter, lf_map, counts, string_length):
    beginning = counts[letter] + lf_map[letter][begin - 1] + 1
    ending = counts[letter] + lf_map[letter][end]
    return(beginning,ending)


def generate_all(input_string, s_array=None, eos="$"):
    letters = set(input_string)
    counts = count_occurences(input_string)

    input_string = "".join([input_string, eos])
    if s_array is None:
        s_array = np.argsort(suffix_array3(input_string))
    bwt = bwt_from_suffix(input_string, s_array)
    lf_map = lf_mapping(bwt)

    for i, j in lf_map.items():
        j.extend([j[-1], 0])
    return letters, bwt, lf_map, counts, s_array


def find_all(string, pattern):
    """
    Parameters
    ----------
    string : TYPE
        DESCRIPTION.
    pattern : TYPE
        DESCRIPTION.

    Yields
    ------
    i : TYPE
        DESCRIPTION.

    """
    i = string.find(pattern)
    while i != -1:
        yield i
        i = string.find(pattern, i+1)
        
def find_v2(search_string, input_string, mismatches_5p=0,mismatches_3p=2, bwt_data=None, s_array=None):
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
    result_new=[]
    #inverted 
    mismatch_flag = 1
    if len(search_string) == 0:
        #return("Empty Query String")
        return []
    if bwt_data is None:
        bwt_data = generate_all(input_string, s_array=s_array)
        
    letters, bwt, lf_map, count, s_array = bwt_data
    
    if len(letters) == 0:
        #return('Empty Query String')
        return []
#If there are some letters not in letters, then disregard immediately
    if not set(search_string) <= letters:
        return []
    length = len(bwt)
    
    
    class Fuzzy(object):
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)
    fuz = [Fuzzy(search_string=search_string, begin=0, end=len(bwt) - 1,
                        mismatches=0)]
    counter = 0
    while len(fuz) > 0:
        p = fuz.pop()
        counter = counter +1
        #List gradually gets smaller and smaller - letter by letter
        searching = p.search_string[:-1]
        #print(searching)
        last = p.search_string[-1]
        #Possible letters the last one could be - mismatching occurs through backtracking remember
        all_letters = [last] if p.mismatches == 0 else letters
        for letter in all_letters:
            
            begin, end = update(p.begin, p.end, letter, lf_map, count, length)
            #mismatch!
            if begin > end:
                mismatch_flag = 0
                #print('we have position ###:', (counter, round(len(search_string)*0.25)))
                if counter <= round(len(search_string)*0.45):
                    #print('True', counter, searching)
                    #print('Mismatch in 3`')
                    search_string_75 = search_string[:round(len(search_string)*0.55)]
                    search_string_25 = search_string[round(len(search_string)*0.55):]
                    
                    #Searching 55% of read - 0 mismatches
                    fuz_5p = [Fuzzy(search_string=search_string_75, begin=0, end=len(bwt) - 1,
                                        mismatches=mismatches_5p)]
                    
                    while len(fuz_5p) > 0:
                        #print('Looking for mismatch in 55% of read')
                        long_read = fuz_5p.pop()
                        
                        search_long = long_read.search_string[:-1]
                       
                        #print(search_long)
                        last_long = long_read.search_string[-1]
                        
                        #Possible letters the last one could be - mismatching occurs through backtracking remember
                        all_letters_long = [last_long] if long_read.mismatches == 0 else letters
                        for letter_long in all_letters_long:
                            
                            begin_long, end_long = update(long_read.begin, long_read.end, letter_long, lf_map, count, length)
                            if begin_long > end_long:
                                #Mismatch in 55% of read found.. exiting')
                                return []
                                
                            if begin_long <= end_long:
                                #Keeps shortening until the string is complete - when this happens, results are extended..
                                if len(search_long) == 0:
                                    #Here is where results get added.
                                    results.extend(s_array[begin_long : end_long + 1])
                                    
                                else:
                                    miss = p.mismatches
                                    if letter_long != last_long:
                                        miss = max(0, long_read.mismatches - 1)
                                    fuz_5p.append(Fuzzy(search_string=search_long, begin=begin_long,
                                                            end=end_long, mismatches=miss))
                    
                    
                    fuz_2 = [Fuzzy(search_string=search_string_25, begin=0, end=len(bwt) - 1,
                                        mismatches=mismatches_3p)]
                    result_new = []
                    while len(fuz_2) > 0:
                        #print('Allowing 2 mismatches in 3`.')
                        second = fuz_2.pop()       
                        searching_2 = second.search_string[:-1]
                        
                        #print(searching_2)
                        last_2 = second.search_string[-1]
                        all_letters_2 = [last_2] if second.mismatches == 0 else letters
                        
                        for base in all_letters_2:
                
                             begin_short, end_short = update(second.begin, second.end, base, lf_map, count, length)
                             if begin_short <= end_short:
                                 
                                 if len(searching_2) == 0:
                
                                     result_new.extend(s_array[begin_short : end_short + 1])
                                     
                                 else:
                                     miss_2 = second.mismatches
                                     if base != last_2:
                                         miss_2 = max(0, second.mismatches - 1)
                                     fuz_2.append(Fuzzy(search_string=searching_2, begin=begin_short,
                                                             end=end_short, mismatches=miss_2))
                           
            if begin <= end:
                next
                if len(searching) != 0:
                    miss = p.mismatches
                    if letter != last:
                        miss = max(0, p.mismatches - 1)
                    
                    fuz.append(Fuzzy(search_string=searching, begin=begin,
                                         end=end, mismatches=miss))

                if len(searching) == 0:
                        results.extend(s_array[begin : end + 1])

    if results == []:
        results = None
    if result_new ==[]:
        result_new = None
        
    return results, result_new, mismatch_flag

class DecimalCounter(Counter):
    """
    A Class and method to calculate and assign counts to multiple
    mapping reads whereby 1 read maps to N loci and each loci receives
    1/N reads.
    """
    def update_division(self, iterable):
        increment = 1 / len(iterable)
        for elem in iterable:
            self[elem] += increment


def process_chunk(input_dict, ref_seq, ids_ref, bwt_data, mismatches_5p=0, mismatches_3p=2):
    """
    Parameters
    ----------
    input_dict : dict
        Dictionary of FASTA.
    ref_seq : str
        Reference sequence to align against.
    ids_ref : list
        Reference sequence IDs.
    bwt_data : list
        BWT created data.

    Returns
    -------
    None.

    """
    #Create empty dictionary for storing alignments
    res_dict = {x:[] for x in input_dict}
    #Run the alignment
    read_counter = 0
    print('Starting alignment..')
    test_dict = DecimalCounter()
    for i in input_dict:
        read_counter +=1
        res_dict[i] = (find_v2(input_dict[i], ref_seq, mismatches_5p=mismatches_5p,mismatches_3p=mismatches_3p, bwt_data=bwt_data))
        factor = 250000
        if read_counter % factor == 0:
            print(f"Aligned {read_counter} sequences...")
            
    print('Processing alignment..')
 
    #Removes reads with no hits
    res_dict = {k: v for k, v in res_dict.items() if any(v)}
    
    #Locates the position of all spaces in the reference
    space_pos = list(find_all(ref_seq, ' '))
    space_pos = list(map(lambda x: x+1, space_pos))
    
    #Refers the hit back to the reference sample
    for key in res_dict.copy():
        res = res_dict[key]
        if (len(res)>0) & (len(res[0]) != 0):
            #Removes reads with 3 mismatches in main part of read.
            if (res[2] == 0) & (len(res[0]) >= 1 ) & (not bool(res[1])):
                res_dict.pop(key)
                continue
        else:
            res_dict.pop(key)
            continue
    
        #Mismatch filter is 0 as mismatch and 1 as no mismatch
        main_pos, pos_3p, mm = res

        #True if there was a mismatch and there is a match in the 3`.
        if bool(pos_3p):
            np_main_pos = np.sort(np.array(main_pos))
            np_pos_3p = np.sort(np.array(pos_3p))

            #Removes any 3` match less than the lowest 5` match - these can't be real.
            np_pos_3p = np_pos_3p[np_pos_3p >= np_main_pos.min()]
            if len(np_pos_3p) == 0:
                res_dict.pop(key)
                continue
            
            #Finds the value immediately greater than each main_pos value
            indices = [bisect_right(np_pos_3p, np_main_pos[x])-1 if bisect_right(np_pos_3p, np_main_pos[x]) >= len(np_pos_3p) else bisect_right(np_pos_3p, np_main_pos[x]) for x in range(len(np_main_pos))]
            main_pos_2 = []
            #For each main_pos value, it takes the closest match at the 3` and looks if there is a space - if so - then it hasn't matched to the same reference seq
            for ind, val in enumerate(indices):
                if (ref_seq[np_main_pos[ind]:np_pos_3p[val]].find(' ') == -1) & (len(ref_seq[np_main_pos[ind]:np_pos_3p[val]]) >0):
                    main_pos_2.append(np_main_pos[ind])
                else:
                    continue
         
            main_pos = main_pos_2
            if len(main_pos) < 1:
                res_dict.pop(key)
                #print('No alignments left, removing ', key)
                continue
        #Multi-mapping reads - will be equally split over number of sites
        if len(main_pos) > 1:
            read_ind = [bisect_left(space_pos, main_reads) for main_reads in main_pos]
            actual_ind_list = []
            #Result changes whether the read aligns exactly to the start of the sequence or elsewhere
            for x in read_ind:
                    if x == len(space_pos):
                        actual_ind = len(space_pos) -1
                        actual_ind_list.append(actual_ind)
#                            print(actual_ind, 'is')
                        continue
                    if space_pos[x] not in main_pos:
                        actual_ind = x-1
                    else:
                        actual_ind=x
 #                       print(actual_ind)
                    actual_ind_list.append(actual_ind)
                
            #Elucidate 5p/3p aligning
            if len(set(actual_ind_list)) != 1:
                ids = [ids_ref[x] for x in read_ind]
                ids = list(set(ids))
                try:
                    subject_seq = [ref_seq[space_pos[actual]:space_pos[actual+1]-1] for actual in actual_ind_list]
                    
                except IndexError:
                    actual_ind_list_new = [x for x in actual_ind_list if x < max(actual_ind_list)]
                    subject_seq_int = [ref_seq[space_pos[actual]:space_pos[actual+1]-1] for actual in actual_ind_list_new]
                    subject_seq_end = [ref_seq[space_pos[actual_ind_list[actual]]:] for actual in range(len(actual_ind_list)-len(actual_ind_list_new))]
                    subject_seq_int.extend(subject_seq_end)
                    subject_seq = subject_seq_int
                names= []
                for i in range(len(ids)):
                    if main_pos[i] + 9 < space_pos[actual_ind_list[i]] + (len(subject_seq[i])/2):
                        orient = '-5p'
                    else:
                        orient = '-3p'
                    int_name = ids[i].split(' MI')[0] + orient
                    names.append(int_name)
                res_dict[key] = names

                try:
                    test_dict.update_division(res_dict[key])
                except NameError:
                    test_dict = DecimalCounter()
                    test_dict.update_division(res_dict[key])
                continue
            else:
                names=[]
                try:
                    space_pos_ind = bisect_left(space_pos, main_pos[0])
                    if space_pos[space_pos_ind] not in main_pos:
                        actual_ind = space_pos_ind-1
                    else:
                        actual_ind = space_pos_ind
                    
                    subject_seq=ref_seq[space_pos[actual_ind]:space_pos[actual_ind+1]-1]
                except:
                    res_dict.pop(key)
                    continue
                names = ids_ref[actual_ind+1].split(' MI')[0] + orient
                res_dict[key] = names
        else:
            try:
                space_pos_ind = bisect_left(space_pos, main_pos[0])
                if space_pos[space_pos_ind] not in main_pos:
                    actual_ind = space_pos_ind-1
                else:
                    actual_ind=space_pos_ind
        
                subject_seq = ref_seq[space_pos[actual_ind]:space_pos[actual_ind+1]-1]

                if main_pos[0] + 9 < space_pos[actual_ind] + (len(subject_seq)/2):
                    orient = '-5p'
                else:
                    orient = '-3p'
            except:
                res_dict.pop(key)
                continue
            #Plus one here because there is no space before the first reference sequence
            names = ids_ref[actual_ind+1].split(' MI')[0] + orient
            res_dict[key] = names
        
            try:
                test_dict.update([res_dict[key]])
        
            except NameError:
                test_dict = DecimalCounter()
                test_dict.update([res_dict[key]])

    
    return test_dict, res_dict
