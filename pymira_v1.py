#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PymiRa - Version 1
Created on Fri Mar 15 22:08:58 2024

@author: zac
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
#from memory_profiler import profile

file_path = sys.argv[1]
ref_path = sys.argv[2]
out_path = sys.argv[3]


def import_fasta(input_path):
    """
    Imports a FASTA file storing the read IDs, the entire FASTA and the 
    indexes of the IDs.
    Parameters
    ----------
    input_path : string
        Path to the FASTA file.

    Returns
    -------
    rna_ids : List
        List of read IDs
    fasta : List
        List of the entire FASTA
    id_indexes : List
        List of the ID indexes
    """
    with open(str(input_path)) as file:
        fasta_raw = file.readlines()
    fasta = [fasta_raw[x].replace('\n','') for x in range(len(fasta_raw))]
    id_indexes = [i for i, s in enumerate(fasta) if '>' in s]
    rna_ids = [fasta[x][1:] for x in id_indexes]
    print(str(input_path) + ' has successfully been uploaded')

    return rna_ids, fasta, id_indexes

def format_fasta_dict(rna_ids,fasta_seq, id_index):
    """
    Creates a dictionary with read IDs as keys and sequences as values
    e.g. {ID:Value}

    Parameters
    ----------
    rna_ids : List
        List of read IDs
    fasta : List
        List of the entire FASTA
    id_indexes : List
        List of the ID indexes

    Returns
    -------
    a_dict : Dict
        A dictionrary with read IDs as keys and sequences as values

    """
    a_dict = {str(rna_ids[x]):
              (str(' '.join(fasta_seq[id_index[x]+1:len(fasta_seq)]).replace(' ','')) 
                               if x == len(id_index)-1 
                               else str(' '.join(fasta_seq[id_index[x]+1:id_index[x+1]]).replace(' ',''))) 
              for x in range(len(id_index))}
    return a_dict

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

def suffix_array(s):

    n = len(s)
    k = 1
    line = to_int_keys_best(s)
    while max(line) < n - 1:
        line = to_int_keys_best(
            [a * (n + 1) + b + 1
             for (a, b) in
             zip_longest(line, islice(line, k, None),
                         fillvalue=-1)])
        k <<= 1
    return line

def inverse_array(l):
    n = len(l)
    ans = [0] * n
    for i in range(n):
        ans[l[i]] = i
    return ans

def bwt_from_suffix(string, s_array=None):
    if s_array is None:
        s_array = suffix_array(string)
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
        s_array = inverse_array(suffix_array(input_string))
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
    To be parallelised across cores for alignment and processing
    - Wrapper for find_v2 function
    
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
    
            #Mismatch filter is 0 as mismatch and 1 as no mismatch
            main_pos, pos_3p, mm = res
    
            #True if there was a mismatch and there is a match in the 3`.
            if bool(pos_3p):
                np_main_pos = np.sort(np.array(main_pos))
                np_pos_3p = np.sort(np.array(pos_3p))
    
                #Removes any 3` match less than the lowest 5` match - these can't be real.
                np_pos_3p = np_pos_3p[np_pos_3p >= np_main_pos.min()]
                if len(np_pos_3p) == 0:
                    continue
                
                #Finds the value immediately greater than each main_pos value
                indices = [bisect_right(np_pos_3p, np_main_pos[x])-1 if bisect_right(np_pos_3p, np_main_pos[x]) >= len(np_pos_3p) else bisect_right(np_pos_3p, np_main_pos[x]) for x in range(len(np_main_pos))]
                main_pos_2 = []
                #For each main_pos value, it takes the closest match at the 3` and looks if there is a space - if so - then it hasn't matched to the same reference seq
                for ind, val in enumerate(indices):
                    if (ref_seq[np_main_pos[ind]:np_pos_3p[val]].find(' ') == -1) & (len(ref_seq[np_main_pos[ind]:np_pos_3p[val]]) >0):
                        main_pos_2.append(np_main_pos[ind])
                    else:
                        break
             
                main_pos = main_pos_2
            
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
                    try:
                        subject_seq = [ref_seq[space_pos[actual]:space_pos[actual+1]-1] for actual in actual_ind_list]
                        
                    except IndexError:
                        actual_ind_list_new = [x for x in actual_ind_list if x < max(actual_ind_list)]
                        subject_seq_int = [ref_seq[space_pos[actual]:space_pos[actual+1]-1] for actual in actual_ind_list_new]
                        subject_seq_end = [ref_seq[space_pos[actual_ind_list[actual]]:] for actual in range(len(actual_ind_list)-len(actual_ind_list_new))]
                        subject_seq_int.extend(subject_seq_end)
                        subject_seq = subject_seq_int
                    names= []
                    for i in range(len(main_pos)):
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

##Main
def main(file_path, ref_path, out_path, mismatches_5p=0, mismatches_3p=2):
    
    st = time.time()
    #Import and format the FASTA
    input_file_ids, input_file_fasta, input_file_indexes = import_fasta(file_path)
    input_file_fasta = list(map(lambda x: x.replace('T', 'U'), input_file_fasta))
    
    #Split up the FASTA into chunks
    div_len = round(len(input_file_ids)/multiprocessing.cpu_count())
    iterators = list(range(0, len(input_file_ids), div_len+1))
    fasta_iterator = list(range(0, len(input_file_fasta), div_len+1))
    last = len(input_file_ids) - iterators[-1]
    counter = 0
    container=[]
    for lo in range(len(iterators)):
        counter += 1
        if lo == len(iterators)-1:
            ref_ids_sub = input_file_ids[iterators[lo]:iterators[lo]+last]
            ref_fasta_sub = input_file_fasta[fasta_iterator[lo*2]:]
            ref_inds_sub = [i for i, s in enumerate(ref_fasta_sub) if '>' in s]
            ref_dict1 = format_fasta_dict(ref_ids_sub, ref_fasta_sub, ref_inds_sub)
            container.append(ref_dict1)
        else: 
            ref_ids_sub = input_file_ids[iterators[lo]:iterators[lo+1]]
            ref_fasta_sub = input_file_fasta[fasta_iterator[lo*2]:fasta_iterator[lo*2+2]]
            ref_inds_sub = [i for i, s in enumerate(ref_fasta_sub) if '>' in s]
            ref_dict1 = format_fasta_dict(ref_ids_sub, ref_fasta_sub, ref_inds_sub)
            container.append(ref_dict1)
    
    #Import and format the reference sequence (miRBase)
    ref_ids, ref_fasta, ref_inds = import_fasta(ref_path)
    ref_fasta = list(map(lambda x: x.replace('T', 'U'), ref_fasta))
    print('Generating reference..')
    ref_dict = format_fasta_dict(ref_ids, ref_fasta, ref_inds)
    ref = [' '.join(list(ref_dict.values()))]
    del ref_dict, ref_fasta, ref_inds
    clean_ref = re.sub(r'[^a-zA-Z\s]', '', str(ref))
    
    #BWT creation
    letters, bwt, lf_map, count, s_array  = generate_all(clean_ref)
    required = [letters,bwt,lf_map,count,s_array]   
    print('Reference creation:',round(time.time()-st,2), '(s)')    
    
    #Alignment
    st = time.time()
    if __name__ == '__main__':
        process_chunk2 = partial(process_chunk, ref_seq = clean_ref, 
                                 ids_ref=ref_ids, bwt_data = required, mismatches_5p=mismatches_5p,
                                 mismatches_3p=mismatches_3p)
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            res_res = pool.map(process_chunk2, container)

    #Results processing - unpack - needs optimising
    emp = []
    for x in range(len(res_res)):
        if len(res_res[x][0]) > 0:
            emp.append(res_res[x][0])
    
    final = emp[0]
    for merge in range(1,len(emp)):
        final.update(emp[merge])
    
    log={}
    for up in range(len(res_res)):
        log.update(res_res[up][1])
    
    #Creating read-alignment log file
    log_file = pd.DataFrame.from_dict(log, orient='index')
    log_file.to_csv(str(out_path + '_' + '_pymira_log.txt'), header=None)
    
    
    #Creating and formatting counts table
    results = pd.DataFrame.from_dict(final, orient='index')
    results.rename(columns={0:'Count'}, inplace=True)
    results = results.sort_values(by=['Count'], ascending=False)
    results.Count = results.Count.astype(int)
    results = results[results.Count>0]
    total = results.Count.sum()
    final_row = pd.DataFrame({'Count':total},index=['TotalCount'])
    results = pd.concat([results,final_row])
    results.index.name = file_path.split('/')[-1]
    print('Alignment and Processing:',time.time()-st)
    results.to_csv(str(out_path + '_pymira_counts.txt'))

main(file_path, ref_path, out_path)
