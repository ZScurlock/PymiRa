# PymiRa - A rapid and accurate classification tool for small non-coding RNAs (sncRNA), including microRNAs
PymiRa utilises a Burrows-Wheeler transformation algorithm to rapidly identify short reads against a reference FASTA file, written in Python.

This can be applied to identify different types of small non-coding RNAs by changing the reference FASTA file sequences are aligned to e.g. for microRNA identification (miRNAs), 
the miRNA hairpin FASTA file derived from miRBase can be used.

# Usage
PymiRa is available for use at this repository but also is available on our webserver found at XXXX.

**Please note the webserver is intended for individual sample use. For larger experiments please install for your own use, thank you.**

## Input
- FASTA file (ideally from small RNA-Sequencing experiment)
- Reference FASTA file to align sncRNAs against. 

## Output
- `.csv` file counts table of the alignments to the Reference FASTA file.
- `.log` file recording the alignment of each read.

By default, PymiRa allows up to two mismatches at the 3' end of a read.


- Counts table of the type and number 
- Use pymira_v1.py
- Best to use in IDE - cmd version in progress
- Default behaviour - main function - mismatch_5p = 0 mismatch_3p = 2.
- Attempts 0 mismatches for entire read. If mismatch then will retry allowing 'mismatch_5p' mismatches at 5' and 'mismatch_3p' at 3'.

`main(fasta_path, reference_path, out_path)`
