### Temporary Notice
This code is provided to journal reviewers and editors for the sole purpose of evaluating the associated manuscript. 
All rights are reserved by the author(s). 
No permission is granted to copy, distribute, modify, or otherwise use this code for any purpose beyond peer review.
 An open-source license (MIT) will be applied upon formal publication of the manuscript, at which point broader use will be permitted.


# PymiRa - A rapid and accurate classification tool for small non-coding RNAs (sncRNA), including microRNAs
PymiRa utilises a Burrows-Wheeler transformation algorithm to rapidly identify short reads against a reference FASTA file, written in Python.

This can be applied to identify different types of small non-coding RNAs by changing the reference FASTA file sequences are aligned to e.g. for microRNA identification (miRNAs), 
the miRNA hairpin FASTA file derived from miRBase can be used.


# Installation
`pip install git+https://github.com/ZScurlock/PymiRa.git`

To ensure PymiRa has been installed correctly, feel free to use the `--self-test` flag. e.g.

```
pymira --self-test

[Importing package] PASS:Import test passed.
[CLI operations] PASS:CLI test successful.
[Alignment] PASS:Alignment test passed
Self-test passed
```

# Usage
PymiRa is available for use at this repository but also is available on our webserver found at www.pymira.co.uk.

**Please note the webserver is intended for individual sample use. For larger experiments please install for your own use, thank you.**


Basic Usage:
`pymira --input_file /path/to/input_file --ref_file /path/to/reference_file --out_path /path/to/results`

## Input
- `FASTA(.gz)` / `FASTQ(.gz)` file (trimmed of any adapters)
- Reference `FASTA(.gz)` / `FASTQ(.gz)` to align sncRNAs against.
- Output basename for filenames e.g. {BASENAME}_pymira_counts.txt


## Output
- `.csv` file counts table of the alignments to the Reference file.
- `log.json` file recording the alignment of each read.
- `alignment_summary.json` file providing a summary of the alignment.

By default, PymiRa allows up to two mismatches at the 3' end of a read. This parameter can be changed using the `--mismatches_3p` flag.

## All options
```
usage: pymira [-h] --input_file INPUT_FILE --ref_file REF_FILE --out_path OUT_PATH [--num_proc NUM_PROC] [--mirna] [--mismatches_3p MISMATCHES_3P] [--self-test]

Align an <input_file> against a <reference_file> allowing up to <mismatches_3p> mismatches at the 3' end of each read. Default is 2.

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE
                        Input sequence file (FASTA / FASTQ, optionally gzipped) to identify small RNAs from.
  --ref_file REF_FILE   Reference sequence file (FASTA / FASTQ, optionally gzipped) used to align the <input_file> against e.g. miRBase hairpin file for microRNA
                        identification
  --out_path OUT_PATH   The dir and basename of the output files e.g. {BASENAME}_pymira_counts.txt
  --num_proc NUM_PROC   Number of processors to be used. Default is 1.
  --mirna               For use with mature miRNA sequences. Restricts alignments to the ends of a precursor (preventing alignment of isomiR / degradation products).
                        Adds '-5p' / '-3p' notation to miRNA alignments. Default is off.
  --mismatches_3p MISMATCHES_3P
                        Number of mismatches allowed in the 3 prime part of a read (Last 45 percent). Default is 2.
  --self-test           Run a quick test to verify installation and correct operation.

```

## Examples
In the `examples/` folder you will find example input, reference and output files.

As an example, we have a file we want to identify miRNAs from.

```
example_input.fa
>read1
UGAGGUAGUAGGUUGUAUAGUU
>read2
CUAUACAAUCUACUGUCUUUC
>read3
CUGUACAGCCUCCUAGCUUUCC
>read4
UGAGGUAGUAGGUUGUGUGGUU
>read5
CUAUACAACCUACUGCCUUCCC
etc..
```

Since we want to identify miRNAs, we will use a reference file such as the miRBase hairpin precursor FASTA file. In our case, our input data is from humans so we have filtered the hairpin reference file to only contain human (hsa) sequences.

```
mirbase_hsa_hairpin.fa
>hsa-mir-3116-2 MI0014129 Homo sapiens miR-3116-2 stem-loop
UAUUGAGUCCCUACUAUGUUCCAGGCACCUACGAUACCCAGUGCCUGGAACAUAGUAGGGACUCAAUA
>hsa-mir-128-1 MI0000447 Homo sapiens miR-128-1 stem-loop
UGAGCUGUUGGAUUCGGGGCCGUAGCACUGUCUGAGAGGUUUACAUUUCUCACAGUGAACCGGUCUCUUUUUCAGCUGCUUC
>hsa-mir-1262 MI0006397 Homo sapiens miR-1262 stem-loop
AUCUACAAUGGUGAUGGGUGAAUUUGUAGAAGGAUGAAAGUCAAAGAAUCCUUCUGGGAACUAAUUUUUGGCCUUCAACAAGAAUUGUGAUAU
>hsa-mir-4767 MI0017408 Homo sapiens miR-4767 stem-loop
ACAUGGGCCCGCGGGCGCUCCUGGCCGCCGCCCGACUUCGGGGCCAGCCGGGGGCAGAGCGCGCGGGAGCCCGAGCGU
>hsa-mir-6829 MI0022674 Homo sapiens miR-6829 stem-loop
CAGCGUGGGCUGCUGAGAAGGGGCAGGGUCCUCCAGCUCAUUCCUCCUGCCUCCUCCGUGGCCUCAG
etc..
```

Now we have everything we need for alignment.

```
pymira --input_file example_input.fa --ref_file mirbase_hsa_hairpin.fa --out_path example_output

Importing input_file: example_input.fa
Importing reference file: mirbase_hsa_hairpin.fa
Generating reference..
Running alignment..
PymiRa Log file written to example_output_pymira_log.json
PymiRa Alignment summary written to example_output_pymira_alignment_summary.json
PymiRa Counts file written to example_output_pymira_counts.txt
Complete

```

Looking at the `example_output_pymira_counts.txt`, we see we get counts at the reference hairpin level.
```
example_input.fa,Count
hsa-mir-520c MI0003158 Homo sapiens miR-520c stem-loop,2
hsa-mir-151a MI0000809 Homo sapiens miR-151a stem-loop,2
hsa-mir-598 MI0003610 Homo sapiens miR-598 stem-loop,2
hsa-mir-10396b MI0033426 Homo sapiens miR-10396b stem-loop,2
```

If want to get mature miRNA counts and to restrict alignment of reads to either end of the reference hairpin, we can use the `--mirna` flag.

```
example_input.fa,Count
hsa-mir-520c-3p,1
hsa-mir-520c-5p,1
hsa-mir-151a-5p,1
hsa-mir-151a-3p,1
hsa-mir-598-5p,1
hsa-mir-598-3p,1
hsa-mir-10396b-5p,1
hsa-mir-10396b-3p,1
```


## Citation





