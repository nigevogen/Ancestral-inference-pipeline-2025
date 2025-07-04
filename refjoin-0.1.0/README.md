# refjoin

A Python package to combines two sequence alignments using a specified shared sequence. 

## Input files and parameters

- `aln_path1`: A relative/absolute path to alignment file 1. 

- `aln_path2`: A relative/absolute path to alignment file 2.

- `ref1_pos`: The order of “reference” sequence in alignment file 1. (0-based index)

- `ref2_pos`: The order of “reference” sequence in alignment file 2. (0-based index)

- `output_path`: A relative/absolute path to an output file. 


## Output file

a FASTA file. The output file includes all sequences from the input files. Sequences in the input files are aligned without alterations in alignments of sequences from an input alignment file. The sequence shared in the input files is duplicated in the output file. 

## Example

Two alignments are given as input files. These alignment files share sequences, `seq_2a` and `seq_2b`, which are identical if gap characters, `-`, are removed. These sequences are used as "reference sequences".  

Integers below nucleotide sequences indicate positions in each alignment so we can refer to particular alignment positions in text. 

### Input: alignment 1

```
>seq_1
A-T
>seq_2a
AG-

123         # Position in input alignment 1
```

### Input: alignment 2

```
>seq_2b
A-G-
>seq_3
ACGC
>seq_4
ACGC

1234        # Position in input alignment 2
```

### Output

```
>seq_1
A--T-
>seq_2a
A-G--
>seq_2b
A-G--
>seq_3
ACG-C
>seq_4
ACG-C

12345       # Position in output
```

The program looks for matched positions by scanning each position of reference sequences in the input files from the left- to right-ends. 

If the current positions indicate the identical nucleotides of reference sequence in both input (e.g., pos 1 of input 1 and input 2), the program considers these positions “matched” and aligns these positions. If both reference sequences indicate nucleotides (non-gap characters) but different nucleotides at the position, this program raises a ValueError. 

If one of the input indicates a nucleotide (e.g., pos 2 of input 1) and the other indicates a gap char at the current position (e.g., pos 2 of input 2), the program considers these positions not matched. This program appends gap chars to an alignment configuration from the input with the reference sequence with gap (e.g., pos 2 of output for the example above). The appended gap characters represent the absence of an nucleotide in any sequence from input 1. The current position stays for input 1 and moves to the next one for input 2. 

If the current positions are aligned to a gap in non-reference sequences (e.g., pos 2 of input 1), this does not affect how the program inserts an output alignment. The program considers the positions matched and combines the alignment configurations from two input. 

If both input alignments contain gaps in reference sequences (e.g., pos 3 of input 1 and pos 4 of input 2), the program constructs two separate alignment configurations (e.g., pos 4 and pos 5 in output). 

## How to use

Python code block below shows how to use this program for the example above. 

```Python
# Import main function from this package.
from refjoin.nucleotide_alignment import refjoin_alignments

# Define values of input arguments.
aln_path1 = 'sample/input_1.aln'
aln_path2 = 'sample/input_2.aln'
ref1_pos = 1
ref2_pos = 0
output_path = 'sample/output.aln'

# Call the main function. 
output_aln_list = refjoin_alignments(aln_path1, aln_path2, ref1_pos, ref2_pos, output_path)
```

Please use `nucleotide_alignment.refjoin_alignments` function for nucleotide sequence alignments and `codon_alignment.refjoin_codon_alignments` for CDS alignments. 

## Dependency

- numpy

## About authors

`refjoin` is written by [Kent Kawashima][KK_github].

<!-- Link below is invisible on  -->
[KK_github]: https://github.com/kentwait 
