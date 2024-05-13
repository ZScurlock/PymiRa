#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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
import sys
import time
file_path = '/media/sf_phd_data/10_mil.fasta'
file_path = '/media/sf_phd_data/drosha_dicer_ko/trimmed/SRR3174960_trimmed.fasta'
ref_path = '/media/sf_phd_data/set2/no_rc_chimira_counts/fixed_counts/human.fa'
def import_fasta(input_path):
    with open(str(input_path)) as file:
        fasta_raw = file.readlines()
    fasta = [fasta_raw[x].replace('\n','') for x in range(len(fasta_raw))]
    id_indexes = [i for i, s in enumerate(fasta) if '>' in s]
    rna_ids = [fasta[x][1:] for x in id_indexes]
    print(str(input_path) + ' has successfully been uploaded')

    return rna_ids, fasta, id_indexes

def format_fasta_dict(rna_ids,fasta_seq, id_index):
    a_dict = {str(rna_ids[x]):(str(' '.join(fasta_seq[id_index[x]+1:len(fasta_seq)]).replace(' ','')) 
                               if x == len(id_index)-1 else str(' '.join(fasta_seq[id_index[x]+1:id_index[x+1]]).replace(' ',''))) 
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
    i = string.find(pattern)
    while i != -1:
        yield i
        i = string.find(pattern, i+1)
        
def find_v2(search_string, input_string, mismatches=0, bwt_data=None, s_array=None):
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
                        mismatches=mismatches)]
    counter = 0
    while len(fuz) > 0:
        p = fuz.pop()
        counter = counter +1
        #if counter ==7:
        #    break
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
                                        mismatches=0)]
                    
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
                                        mismatches=2)]
                    result_new = []
                    while len(fuz_2) > 0:
                        #print('Allowing 2 mismatches in 3`.')
                        second = fuz_2.pop()
                #         
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
    def update_division(self, iterable):
        increment = 1 / len(iterable)
        for elem in iterable:
            self[elem] += increment
            
    
##Main
st = time.time()
input_file_ids, input_file_fasta, input_file_indexes = import_fasta(file_path)
input_file_fasta = list(map(lambda x: x.replace('T', 'U'), input_file_fasta))
#input_dict = format_fasta_dict(input_file_ids, input_file_fasta, input_file_indexes)
#div_len = len(input_file_ids)/multiprocessing.cpu_count()
iterators = list(range(0, len(input_file_ids), 702917))
fasta_iterator = list(range(0, len(input_file_fasta), 702917))
last = len(input_file_ids) - iterators[-1]
fasta_last = len(input_file_fasta) - fasta_iterator[-1]
counter = 0
container=[]
for lo in range(len(iterators)):
    counter += 1
    if lo == len(iterators)-1:
        ref_ids_sub = input_file_ids[iterators[lo]:iterators[lo]+last]
        #ref_fasta_sub = ref_fasta[fasta_iterator[lo*2]:fasta_iterator[lo*2+1]+fasta_last]
        ref_fasta_sub = input_file_fasta[fasta_iterator[lo*2]:]
        ref_inds_sub = [i for i, s in enumerate(ref_fasta_sub) if '>' in s]
        ref_dict1 = format_fasta_dict(ref_ids_sub, ref_fasta_sub, ref_inds_sub)
        
        #ref1 = [' '.join(list(ref_dict1.values()))]
        container.append(ref_dict1)
        # clean_ref1 = re.sub(r'[^a-zA-Z\s]', '', str(ref1))
        # print('Generating the bwt number', counter)
        # letters, bwt, lf_map, count, s_array  = generate_all(clean_ref1)
        # required = [clean_ref1,letters,bwt,lf_map,count,s_array]

        # with open('hsa' + '_' + str(counter) + '_ref_pickle2.pkl', 'wb') as file:
        #      pickle.dump(required, file)
        # print('hsa' + '_' + str(counter) + '_ref_pickle2.pkl saved.')
    else: 
        ref_ids_sub = input_file_ids[iterators[lo]:iterators[lo+1]]
        ref_fasta_sub = input_file_fasta[fasta_iterator[lo*2]:fasta_iterator[lo*2+2]]
        ref_inds_sub = [i for i, s in enumerate(ref_fasta_sub) if '>' in s]
        ref_dict1 = format_fasta_dict(ref_ids_sub, ref_fasta_sub, ref_inds_sub)
   
        #ref1 = [' '.join(list(ref_dict1.values()))]
        container.append(ref_dict1)
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

print(time.time()-st)    
def process_chunk(input_dict):
    """
    Process a chunk of the FASTA file.
    Returns a Counter object representing the frequencies of each sequence.
    """
    #Import human.fa reference
       
    #Create empty dictionary for each read
    res_dict = {x:[] for x in input_dict}
    #Run the alignment
    read_counter = 0
    print('Starting alignment..')
    for i in input_dict:
        read_counter +=1
        res_dict[i] = (find_v2(input_dict[i], clean_ref, mismatches = 0, bwt_data=required))
        factor = 250000
        if read_counter % factor == 0:
            print(f"Aligned {read_counter} sequences...")
            
    print('Processing alignment..')
    #b = time.time()-uptoalign
    #print('Alignment time:',b-a) 
    #Removes reads with no hits
    res_dict = {k: v for k, v in res_dict.items() if any(v)}
    
    #Locates the position of all spaces in the reference
    space_pos = list(find_all(clean_ref, ' '))
    space_pos = list(map(lambda x: x+1, space_pos))
    
    #Refers the hit back to the reference sample
    count=0
    for key in res_dict.copy():
        count +=1
        factor = 1000
        if count % factor ==0:
            print(f"PROCESSED {count} sequences...")
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
                    if (clean_ref[np_main_pos[ind]:np_pos_3p[val]].find(' ') == -1) & (len(clean_ref[np_main_pos[ind]:np_pos_3p[val]]) >0):
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
                    ids = [ref_ids[x] for x in read_ind]
                    try:
                        subject_seq = [clean_ref[space_pos[actual]:space_pos[actual+1]-1] for actual in actual_ind_list]
                        
                    except IndexError:
                        actual_ind_list_new = [x for x in actual_ind_list if x < max(actual_ind_list)]
                        subject_seq_int = [clean_ref[space_pos[actual]:space_pos[actual+1]-1] for actual in actual_ind_list_new]
                        subject_seq_end = [clean_ref[space_pos[actual_ind_list[actual]]:] for actual in range(len(actual_ind_list)-len(actual_ind_list_new))]
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
            
                    subject_seq = clean_ref[space_pos[actual_ind]:space_pos[actual_ind+1]-1]
    
                    if main_pos[0] + 9 < space_pos[actual_ind] + (len(subject_seq)/2):
                        orient = '-5p'
                    else:
                        orient = '-3p'
                except:
                    res_dict.pop(key)
                    continue
                #Plus one here because there is no space before the first reference sequence
                #print('REF:',ref_dict[ref_ids[actual_ind+1]])
                #print('QUY:',subject_seq)
                names = ref_ids[actual_ind+1].split(' MI')[0] + orient
                res_dict[key] = names
            
                try:
                    test_dict.update([res_dict[key]])
            
                except NameError:
            #                    final_dict = res_dict
                #test_dict = Counter(res_dict.values())
                    test_dict = DecimalCounter()
                    test_dict.update([res_dict[key]])
    #print('Processing time:', time.time()-uptoalign-b)
    #print(file_path,'Overall:', time.time()-uptoalign)
   
    #Format results
    return test_dict
    """
    Ignore this section for parallel
    
    
    results = pd.DataFrame.from_dict(test_dict, orient = 'index')
    results.index.name = file_path.split('/')[-1]
    results.rename(columns={0:'Count'}, inplace=True)
    #results.Count = [('%f' % results.Count.iloc[x]).rstrip('0').strip('.') for x in range(len(results))]
    results.Count = [f'{i:g}' for i in round(results['Count'].sort_values(ascending=False))]
    results.Count = results.Count.astype(int)
    results = results[results.Count>0]
    total = results.Count.sum()
    final_row = pd.DataFrame({'Count':total},index=['TotalCount'])
    results = pd.concat([results,final_row])
    #print('Completed for ', file_path)
    #results.to_csv(out_path)
    return results

    """
st = time.time()
#res_res = []
num_processors = 6
with multiprocessing.Pool(processes=num_processors) as pool:
    res_res = pool.map(process_chunk, container)

merged_dict = Counter()
final = res_res[0]
for merge in range(1,len(res_res)):
    final.update(res_res[merge])



results = pd.DataFrame.from_dict(final, orient='index')
results.index.name = file_path.split('/')[-1]
results.rename(columns={0:'Count'}, inplace=True)
#results.Count = [('%f' % results.Count.iloc[x]).rstrip('0').strip('.') for x in range(len(results))]
results = results.sort_values(by=['Count'], ascending=False)
#results.Count = [f'{i:g}' for i in round(results.sort_values(by=['Count'],ascending=False))]
#results.Count = [f'{i:g}' for float(i) in round(results)]
results.Count = results.Count.astype(int)
results = results[results.Count>0]
total = results.Count.sum()
final_row = pd.DataFrame({'Count':total},index=['TotalCount'])
results = pd.concat([results,final_row])
print(time.time()-st)
#print('Completed for ', file_path)
#results.to_csv(out_path)
results.to_csv('/media/sf_phd_data/drosha_dicer_ko/trimmed/SRR3174960_trimmed.fasta_counts.txt')
