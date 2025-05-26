# p123

A pipeline to infer lineage-specific codon changes using BTW output.  
Step1, step2 and step3 are codes to infer fixations and polymorphims.  
Please check readme files in each directory for the detail.  
This pipeline is developed by Tomotaka Matsumoto. 

This pipeline requires `codonpaths` Python package which is included in ___step3__run_codon_path_and_infer_fixations_polymorphisms/codonpaths-0.3.0.  
`codonpaths` generates all possible paths to change from one codon to another codon that have different nucleotides at multiple positions (e.g., AAA -> AGG). Each step of a path is a nucleotide change at only one position (e.g., AAA -> AGA). Only paths with minimal numbers of steps are considered (i.e., multiple hits at a position are not considered). This package also calculates weighted probabilities of paths using the number of synonymous and nonsynonymous changes in a path and a given ratio of synonymous-to-nonsynonymous substitution rates. This package is developed by Kent Kawashima. 
