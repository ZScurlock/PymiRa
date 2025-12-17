import argparse
import PymiRa.functions as pym
import sys
import re
from functools import partial
import multiprocessing
from collections import defaultdict

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument('--self-test', action="store_true")    
    pre_args,remaining = pre_parser.parse_known_args(argv)
    
    if pre_args.self_test:
        from PymiRa.self_test import test_pymira
        ok, results = test_pymira.run_self_test()
        for name, passed, msg in results:
            print(f"[{name}] {'PASS' if passed else 'FAIL'}:{msg}")
        if ok:
            print('Self-test passed')
            return 0
        else:
            print('Self-test FAILED')
            return 1
    
    
    parser = argparse.ArgumentParser(
        description="Align an <input_file> against a <reference_file> allowing up to <mismatches_3p> mismatches at the 3' end of each read. Default is 2. Outputs sequence counts."
    )
    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Input sequence file (FASTA / FASTQ, optionally gzipped) to identify small RNAs from.",
    )
    parser.add_argument(
        "--ref_file",
        type=str,
        required=True,
        help="Reference sequence file (FASTA / FASTQ, optionally gzipped) used to align the <input_file> against e.g. miRBase hairpin file for microRNA identification",
    )
    parser.add_argument(
        "--out_path",
        type=str,
        required=True,
        help="The dir and basename of the output files e.g. {BASENAME}_pymira_counts.txt",
    )
    parser.add_argument(
        "--num_proc",
        type=int,
        default=1,
        help="Number of processors to be used. Default is 1.",
    )
    
    parser.add_argument(
        '--mirna',
        action='store_true',
        default=False,
        help="""For use with mature miRNA sequences.\n Restricts alignments to the ends of a precursor (preventing alignment of isomiR / degradation products).
        Adds '-5p' / '-3p' notation to miRNA alignments.Default is off.
        """
        )
   
    parser.add_argument(
        "--mismatches_3p",
        type=int,
        default=2,
        help="Number of mismatches allowed in the 3 prime part of a read (Last 45 percent). Default is 2.",
    )
    parser.add_argument(
        "--self-test", action="store_true",default=False,
        help="Run a quick test to verify installation and correct operation")
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    
        
    # Import and format the FASTA
    print(f"Importing input_file: {args.input_file}")
    input_file_dict, input_file_ids = pym.upload_file(args.input_file)
    
    seq_dict = defaultdict(list)
    for name,seq in input_file_dict.items():
        seq_dict[seq].append(name)
    
    int_dict = {k:len(v) for k,v in seq_dict.items()}
    
    container = pym.split_fasta_dict(int_dict, int(args.num_proc))

    # Import and format the reference sequence (miRBase)
    print(f"Importing reference file: {args.ref_file}")
    ref_dict, ref_ids = pym.upload_file(args.ref_file)

    print("Generating reference..")
    ref = [" ".join(list(ref_dict.values()))]
    clean_ref = re.sub(r"[^a-zA-Z\s]", "", str(ref))
    clean_ref += ' '
    clean_ref = ' ' + clean_ref
    
    # BWT creation
    letters, bwt, lf_map, count, s_array = pym.generate_all(clean_ref)
    required = [letters, bwt, lf_map, count, s_array]

    # Alignment
    print("Running alignment..")
    process_chunk2 = partial(
        pym.process_chunk,
        ref_seq = clean_ref,
        ids_ref = ref_ids,
        bwt_data = required,
        integer_dict = int_dict,
        mismatches_3p = args.mismatches_3p,
        mirna_flag = args.mirna
    )
    with multiprocessing.Pool(processes=args.num_proc) as pool:
        res_res = pool.map(process_chunk2, container)

    # Process result
    try:
        final = pym.merge_results(res_res)
        pym.create_log(res_res,seq_dict,args.out_path)
        pym.generate_results_files(final, args.input_file, args.ref_file, args.out_path,
                                   args.mirna, args.mismatches_3p,input_file_dict)
        print('Complete')
    except IndexError:
        print('No read sequences were aligned in <input_file> so there are no results.')
        
if __name__ == "__main__":
    main()
