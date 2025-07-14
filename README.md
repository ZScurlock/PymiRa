# PymiRa - A rapid and accurate classification tool for small non-coding RNAs (sncRNA), including microRNAs
PymiRa utilises a Burrows-Wheeler transformation algorithm to rapidly identify short reads against a reference FASTA file, written in Python.

This can be applied to identify different types of small non-coding RNAs by changing the reference FASTA file sequences are aligned to e.g. for microRNA identification (miRNAs), 
the miRNA hairpin FASTA file derived from miRBase can be used.

#Installation
`pip install git+https://github.com/ZScurlock/PymiRa.git`

# Usage
PymiRa is available for use at this repository but also is available on our webserver found at XXXX.

`pymira --input_fasta /path/to/input_fasta --ref_fasta /path/to/ref_fasta --out_path /path/to/results`

**Please note the webserver is intended for individual sample use. For larger experiments please install for your own use, thank you.**

## Input
- FASTA file (ideally from small RNA-Sequencing experiment)
- Reference FASTA file to align sncRNAs against.
- Output file names


## Output
- `.csv` file counts table of the alignments to the Reference FASTA file.
- `.log` file recording the alignment of each read.

By default, PymiRa allows up to two mismatches at the 3' end of a read. If you wish to disable this feature, please set the `mismatches_3p` to `0` at your own discretion.

#All options
usage: pymira [-h] --input_fasta INPUT_FASTA --ref_fasta REF_FASTA --out_path OUT_PATH [--num_proc NUM_PROC] [--mismatches_5p MISMATCHES_5P] [--mismatches_3p MISMATCHES_3P]

Align a <input_fasta> against a <reference_fasta> file, to obtain a counts and log file.

optional arguments:
  -h, --help            show this help message and exit
  --input_fasta INPUT_FASTA
                        An input FASTA file to identify small RNAs from.
  --ref_fasta REF_FASTA
                        A reference FASTA file used to align the <input_fasta> against e.g. miRBase hairpin file for microRNA identification
  --out_path OUT_PATH   The dir and basename of the output files
  --num_proc NUM_PROC   Number of processors to be used. Default is 4.
  --mismatches_5p MISMATCHES_5P
                        Number of mismatches allowed in the 5 prime part of a read (First 55 percent). Default is 0, change at your own discretion
  --mismatches_3p MISMATCHES_3P
                        Number of mismatches allowed in the 3 prime part of a read (Last 45 percent). Default is 2.

