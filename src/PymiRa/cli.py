import argparse
import PymiRa.functions as pym
import sys
import re
from functools import partial
import multiprocessing
import pandas as pd

def main():
    parser=argparse.ArgumentParser(description='Align an <input_fasta> against a <reference_fasta> file, to obtain a counts and log file.')
    parser.add_argument('--input_fasta',type=str,required=True,help='An input FASTA file to identify small RNAs from.')
    parser.add_argument('--ref_fasta',type=str,required=True,help='A reference FASTA file used to align the <input_fasta> against e.g. miRBase hairpin file for microRNA identification')
    parser.add_argument('--out_path', type=str, required=True,help='The dir and basename of the output files')
    parser.add_argument('--num_proc', type=int, default=4,help='Number of processors to be used. Default is 4.')

    parser.add_argument('--mismatches_5p', type=int, default=0, help='Number of mismatches allowed in the 5 prime part of a read (First 55 percent). Default is 0, change at your own discretion')
    parser.add_argument('--mismatches_3p',type=int,default=2, help='Number of mismatches allowed in the 3 prime part of a read (Last 45 percent). Default is 2.')
    

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args=parser.parse_args()
    #if len(sys.argv) == 1:
    #    parser.print_help()
    
    print(f'{args.input_fasta}')
    


    num_proc = int(args.num_proc)
#    st = time.time()
    #Import and format the FASTA
    input_file_dict, input_file_ids = pym.parse_fasta(args.input_fasta)
    #print(input_file_dict[list(input_file_dict.keys())[0]])
    #input_file_fasta = list(map(lambda x: x.replace('T', 'U'), input_file_fasta))
    container = pym.split_fasta_dict(input_file_dict, int(args.num_proc))

    #Import and format the reference sequence (miRBase)
    ref_dict, ref_ids  = pym.parse_fasta(args.ref_fasta)
    #ref_fasta = list(map(lambda x: x.replace('T', 'U'), ref_fasta))
    print('Generating reference..')
    #ref_dict = pym.format_fasta_dict(ref_ids, ref_fasta, ref_inds)
    ref = [' '.join(list(ref_dict.values()))]
    del ref_dict
    clean_ref = re.sub(r'[^a-zA-Z\s]', '', str(ref))

    #BWT creation
    letters, bwt, lf_map, count, s_array  = pym.generate_all(clean_ref)
    required = [letters,bwt,lf_map,count,s_array]
   # print('Reference creation:',round(time.time()-st,2), '(s)')

    #Alignment
 #   st = time.time()
    process_chunk2 = partial(pym.process_chunk, ref_seq = clean_ref, ids_ref=ref_ids, bwt_data = required, mismatches_5p=args.mismatches_5p, mismatches_3p=args.mismatches_3p)
    with multiprocessing.Pool(processes=args.num_proc) as pool:
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
    log_file.to_csv(str(args.out_path + '_pymira_log.txt'), header=None)


    #Creating and formatting counts table
    results = pd.DataFrame.from_dict(final, orient='index')
    results.rename(columns={0:'Count'}, inplace=True)
    results = results.sort_values(by=['Count'], ascending=False)
    results.Count = results.Count.astype(int)
    results = results[results.Count>0]
    total = results.Count.sum()
    final_row = pd.DataFrame({'Count':total},index=['TotalCount'])
    results = pd.concat([results,final_row])
    results.index.name = args.input_fasta.split('/')[-1]
  #  print('Alignment and Processing:',time.time()-st)
    results.to_csv(str(args.out_path + '_pymira_counts.txt'))

if __name__ == "__main__":
    main()
