# baseml-btw-change_counts

This pipeline estimates ancestral states and counts of nucleotide changes for mutation types. Our study (Yamashita et al 2025) used these codes to study Drosophila but they are applicable to any data set of aligned sequences from within- and between-species. An assumed gene tree is required and the parameter-rich model necessitates large data sets.

## Process overview

### 1. Clean input alignment files
- exclude sites/codons aligned to gaps, N chars, and other specified sites/codons. 
- exclude sequences that are not employed for ancestral inference. 

### 2. Run BASEML (included in PAML, Yang 2007)
- make “collapse” sequences. (See Fig. 2 of Matsumoto and Akashi 2018)
- execute BASEML. ([baseml installation][PAML_github] is required before starting this pipeline) 
- See [PAML manual and wiki in Github repository][PAML_github] for detail. 

### 3. Run BTW
- adjust probabilities of ancestral states using theoretically expected SFS for neutral mutations.
- construct SFS using ancestral states and the adjusted probabilities. 
- use the estimated SFS for adjusting the probabilities of ancestral states. Repeat the process of SFS construction and probability adjustment for a given number of iterations. 
- See Matsumoto and Akashi (2018) about the BTW method. The details of computational processes and output files are described in "btw/BTW_instruction_TM_ver1.1.pdf". 

### 4. Infer nucleotide change counts
- Compute the sum of probabilities across internal node nucleotide configurations showing different nucleotide states at ancestral and derived nodes. This is an implementation of AWP method introduced in Matsumoto and Akashi (2018). 
- Output summary: total # of sites and # of variable sites in the clean alignment, base composition at each node, fixation counts and SFS (counts) for each mutation type, and terminal node nucleotide configurations and their counts in the clean alignment. 
- Perform bootstrap resampling of estimated nucleotide change counts in units of input alignment files (often CDS or introns) and pool counts across sites for each branch for each replicate. 


## Input files

- control file: “0.ctl.ini”
The template file is in this folder. 
- a folder of alignment FASTA files. 
The folder path should be specified in the control file as “input_aln_dir” parameter.
“Marker sequences” can be included to indicate sites to exclude (see [below](#marker_seq_bookmark)). 
- a file for a list of alignment file names from the input alignment folder. 
The file path should be specified in the control file as “input_aln_list_file” parameter. 
e.g., “0.ms_auto.ms_short.cutLowrec.itemlist” in intron alignment folder (“m14s21ye_r624_intron_aln_201218HYc”) in [Zenodo][zenodo_link]. 
- a newick file to indicate an assumed species tree. 
This file is passed to BASEML and should follow the format required by BASEML (e.g., how to specify branches to share transition parameters). 

Below is from p.15 of PAML manual (version 4.9j, included in btw folder):  
>_Branch or node labels._ Some models implemented in baseml and codeml allow several groups of branches on the tree, which are assigned different parameters of interest. For example, in the local clock models (clock = 2 or 3) in baseml or codeml, you can have, say, 3 branch rate groups, with low, medium, and high rates respectively. Also the branch-specific codon models (model = 2 or 3 for codonml) allow different branch groups to have different ωs, leading to so called “two- ratios” and “three-ratios” models. All those models require branches or nodes in the tree to be labeled. Branch labels are specified in the same way as branch lengths except that the symbol “#” is used rather than “:”. The branch labels are consecutive integers starting from 0, which is the default and does not have to be specified. For example, the following tree  
>`((Hsa_Human, Hla_gibbon) #1, ((Cgu/Can_colobus, Pne_langur), Mmu_rhesus), (Ssc_squirrelM, Cja_marmoset));`  
>is from the tree file examples/lysozyme/lysozyme.trees, with a branch label for fitting models of different ω ratios for branches. The internal branch ancestral to human and gibbon has the ratio ω1, while all other branches (with the default label #0) have the background ratio ω0. This fits the model in table 1C for the small data set of lysozyme genes in Yang (1998). See the readme file in the examples/lysozyme/ folder.


## Output files

This pipeline will create three folders under a specified output folder: “1_alignment_processing”, “2_BTW” and “3_SFS”. 

1_alignment_processing/  
0.alignment.ini  
0.read_ctl.txt  
1_unfiltered_with_markers.zip  
2_filtered_without_markers.zip  
2_filtered_without_markers.csv  

2_BTW/  
run_BTW_pipeline.sh  
____codes  
_AI_for_collapse	# copy from this folder  
_collapse_seqs_folder  # “collapse-pair sequences” in Fig. 2 of Matsumoto and Akashi 2018  
_concatenated_original_alignment_folder  
_seqs_folder  
_seq_list_to_be_concatenated  
_tree_folder/  
sample.trees		# given newick file  
BASEML_result/  
baseml_output/	 	# output folder of BASEML  
iterative_BTWest/  
anc_site_probs_all_node.txt	# ancestral reconstruction for each alignment position  

3_SFS/  
anc_state_summary.txt	# input site counts, base composition at each node, SFS and fixation counts  
joint_probs.txt.run_info.txt  
joint_probs.txt		# nucleotide change counts for each branch for each alignment position  

## How to run the pipeline

### Preparation

Change `package_dir` value to the location of this folder on your computer. The variable is defined in line 14 of `run.sh` and in line 8 of `utils/make_output_dirs.sh` like below,

```
package_dir=/Users/<username>/<thisfolder>
```

### Run

1. Prepare a control file. A template file, “0.ctl.ini”, is available in this repository.

2. Create an output folder (specified in the control file) if it does not exists. 

3. In Terminal, go to baseml-btw-change_count folder on your computer.

4. Run this command

```
run.sh <control file>
```

## Parameters in the control file

- `site_type`: codon or nucleotide. If codon is given, input sequence length has to be multiple of 3. Only the codon 3rd positions are used for inferring ancestral states and nucleotide changes.
- `sample_order`: The order of codon/nucleotide sequences. This parameter is required to guarantee the consistency of the species order across alignment files. Sequence names are separated by “,”. Regular expression can be used to specify sequence names if sequence names are different among alignment files for a given species (e.g., “Dyak_\d+” and “Dere_\d+” below). 
    - Example from Dmel subgroup dataset (zenodo link):  
    RG18N$, RG2$, RG19$, RG22$, RG24$, RG25$, RG28$, RG32N$, RG34$, RG36$, RG38N$, RG3$, RG5$, RG9$, MD03$, MD06$, MD105$, MD106$, MD146$, MD15$, MD197$, MD199$, MD201$, MD221$, MD224$, MD225$, MD233$, MD235$, MD238$, MD243$, MD251$, MD255$, MD63$, MD72$, MD73$, Dyak_\d+, Dere_\d+
- `marker_kw` (optional): A suffix string of the names of marker sequences. (default: marker)
- `markers_to_use` (optional): The names of marker sequences to use for alignment filtering. The names have to separated by “, ”. 
    - <a name="marker_seq_bookmark">About marker sequence</a>  
    A marker sequence can contain “N” and “C” characters. Sites aligned to “N” are excluded whereas sites aligned to “C” are retained in the cleaned alignment files which are input of BASEML. 
- `genetic_code_type`: a or b. In a-type genetic code, all six fold amino acids are in one class, whereas in b-type Ser codons are separated to four-fold (TCN) and two-fold (AGY) degenerate classes.
- `only_conserved_aa`: True or False. Whether to retain only “conserved” codon positions in the cleaned alignment files. Here “conserved” is defined as identical amino acids at all terminal nodes. 
- `mt2_seg_allele_filter_mode`: third_pos or any.  In the BTW process, within-species SNPsvariations are represented as two sequences (“collapse-pair sequences” in Fig. 2 of Matsumoto and Akashi 2018). Therefore, aligned positions at which >2where 3 or 4 types of nucleotides are segregating within in thea population sample need to be excluded. This option controls which codon positions the code should check for the number of nucleotide types segregating in a population.   
If the pipeline outputs are used to reconstruct codon changes (including nonsynonymous changes), codon 1st, 2nd and 3rd positions are employed to the input of this pipeline. For such analysis, this filtering step needs to check the number of within-species nucleotide types at each codon position and specify “any” for this option in that case. If the analysis is only for codon 3rd positions, use “third_pos”. 

- `trim_kw`: None or <seq_name>:<5p_num>..<3p_num>  
    - where  
    seq_name is a name of sequence to use for site counting  
    5p_num is # of sites to exclude from 5′ end  
    3p_num is # of sites to exclude from 3′ end  
    - If None is given, no trimming is conducted (default: None).  
    Example to specify trimming:
    RG18N:10..30  
    for intron analyses in Yamashita et al 2025 (excluding 10 bases from 5′ end and 30 bases from 3′ end). RG18N (within-species line sequence extracted using Dmel reference annotation) is used for counting bases. 

- `poly_info`: Information for within-species lines about how sequences in input alignment files correspond to the assumed species tree and collapse sequences in the cleaned alignments.
    - Six fields separated by tab (\t):  
    1st: a prefix string of sample sequences  
    2nd: a prefix string of collapse sequences  
    3rd: the order of the collapse sequences that appear in the clean alignments (separated by ",")  
    4th: the index of ancestral node assigned by BASEML  
    5th the order of the first sequence from the population in the input alignment (0-based index).  
    6th: the number of sequences from a population
    
    - Example from the Dmel subgroup dataset:  
    `RG	RG_collapse	3,4	9	0	14  `  
	`MD	MD_collapse	5,6	10	14	21  `

- `single_allele_sp`: Correspondence between sequence order in the input alignment and sequence order in the cleaned alignment for species without within-species lines (only single representatives are included in the input alignment).  
For a given species, the order in the input alignment and the order in the clean alignment are separated by “:”. A part for another species may follow with a separator “, ”. 
    - Example from the Dmel subgroup dataset:  
    `35:1, 36:2`
    - The 35th sequence of input alignment (Dyak sequence) will be placed at the first (top) in the clean alignment. Similar information is listed for Dere. Dere sequence is placed at the 36th in input alignment and at the 2nd in the clean alignment. 

- `mutation_category_num`: The number of nucleotide change categories. (default: 12). Currently this parameter does control an option for pooling because this control file does not modify nucleotide change categories defined in BTW perl code. We may change this parameter to accept user-defined pooling in the future. For now please use 12. 

- `iteration_num`: The number of iterations for BTW. (default: 5) 

- `internal_node_num`: The number of internal nodes of the given species tree. 

- `assert_file_exist`: True or False (default: True)

- `verbose_output`: True or False (default: True)

Parameters not listed here are not incorporated for controlling the process yet. Those parameters are currently defined in shell scripts. 

## Dependencies

- baseml (see [PAML Github repository][PAML_github])
- alignmentrs (see [alignmentrs Github repository][alignmentrs_github])
- pandas
- numpy

## Contributions

- [Ziheng Yang][ZY_github]
Developed the BASEML option to specify branches to share transition parameters.
- [Tomotaka Matsumoto][TM_github]
Developed a pipeline starting from clean alignment and executing BASEML and BTW. 
- [Kent Kawashima][KK_github]
Developed `alignmentrs` which made the alignment filtering process faster. 
- [Haruka Yamashita][HY_github]
Developed Python functions that uses Alignment data structure from `alignmentrs` for alignment filtering. 
Developed Python functions to count nucleotide changes. 
Organized the processes from cleaning input alignment files to inferring nucleotide changes into a pipeline and wrote this README file. 

## References

Yamashita, H., Matsumoto, T., Kawashima, K., Abdulla Daanaa, H. S., Yang, Z., & Akashi, H. (2025). Dinucleotide preferences underlie apparent codon preference reversals in the Drosophila melanogaster lineage. Proceedings of the National Academy of Sciences, 122(21), e2419696122. https://doi.org/10.1073/pnas.2419696122

Yang, Z. (2007). PAML 4: Phylogenetic analysis by maximum likelihood. Molecular Biology and Evolution, 24(8), 1586–1591. https://doi.org/10.1093/molbev/msm088

Matsumoto, T., & Akashi, H. (2018). Distinguishing Among Evolutionary Forces Acting on Genome-Wide Base Composition: Computer Simulation Analysis of Approximate Methods for Inferring Site Frequency Spectra of Derived Mutations in Recombining Regions. G3, 8(5), 1755–1769. https://doi.org/10.1534/g3.117.300512

<!-- Links below do not appear. -->
[PAML_github]: https://github.com/abacus-gene/paml 
[zenodo_link]: https://doi.org/10.5281/zenodo.15274324
[alignmentrs_github]: https://github.com/kentwait/alignmentrs
[ZY_github]: https://github.com/ziheng-yang 
[TM_github]: https://github.com/tomotakamatsumoto 
[KK_github]: https://github.com/kentwait 
[HY_github]: https://github.com/yamasampo

