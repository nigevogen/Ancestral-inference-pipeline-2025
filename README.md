# Ancestor-inference-pipeline-2025

A collection of codes used to infer lineage-specific nucleotide/codon changes for polymorphism/divergence data. These codes can be used to reproduce results in Yamashita et al. 2025 PNAS in combination with sequence alignment data (available at [Zenodo][1]). 

## Packages

This repository consists of three packages: `refjoin`, `ppda`, `codonpaths`, and `p123`. 

### refjoin

Joins sequence alignments that are generated separately and contain an identical sequence. This package is developed by Kent Kawashima. `numpy` is required to use this package. 

### ppda

Infer internal node nucleotide configurations on an assumed species tree for given sequence alignments. BTW is developed by Tomotaka Matsumoto. Other Python functions (e.g., for filtering aligned sites) and shell scripts to automate processes are developed by Haruka Yamashita. 

### codonpaths

Generate all possible paths from one codon to another codon that have different nucleotides at multiple positions (e.g., AAA -> AGG). Each step of a path is a nucleotide change at only one position (e.g., AAA -> AGA). Only paths with minimal numbers of steps are considered (i.e., multiple hits at a position is not considered). This package also calculate weighted probabilities of paths using the number of synonymous and nonsynonymous changes in a path and a given ratio of synonymous-to-nonsynonymous substitution rates. 

## References

Yamashita, H., Matsumoto, T., Kawashima, K., Abdulla Daanaa, H. S., Yang, Z., & Akashi, H. (2025). Dinucleotide preferences underlie apparent codon preference reversals in the Drosophila melanogaster lineage. Proceedings of the National Academy of Sciences, 122(21), e2419696122. https://doi.org/10.1073/pnas.2419696122

Yamashita, H., Matsumoto, T., Kawashima, K., Daanaa, H. S. A., Yang, Z., & Akashi, H. (2024). Recent codon preference reversals in the Drosophila melanogaster lineage (p. 2024.10.10.617326). bioRxiv. https://doi.org/10.1101/2024.10.10.617326

<!-- Links: below are not visible -->
[1]: https://doi.org/10.5281/zenodo.15274324

