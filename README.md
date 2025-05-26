# Ancestor-inference-pipeline-2025

A collection of codes used to infer lineage-specific nucleotide/codon changes for polymorphism/divergence data. These codes can be used to reproduce results in Yamashita et al. 2025 PNAS in combination with sequence alignment data (available at [Zenodo][1]). 

## Packages

This repository consists of one Python package, `refjoin`, and two analysis pipelines, `ppda`, and `p123`. 

### refjoin

Joins sequence alignments that are generated separately and contain an identical sequence. This package is developed by Kent Kawashima. `numpy` is required to use this package. 

### ppda

Infer internal node nucleotide configurations on an assumed species tree for given sequence alignments. BTW is developed by Tomotaka Matsumoto. See Matsumoto and Akashi 2018 (G3). Other Python functions (e.g., for filtering aligned sites) and shell scripts to automate processes are developed by Haruka Yamashita. 

### p123

A pipeline to infer lineage-specific codon changes using BTW output. The pipeline is mainly developed by Tomotaka Matsumoto and includes a Python package `codonpath` developed by Kent Kawashima. Haruka Yamashita contributed to debugging. 

## References

Yamashita, H., Matsumoto, T., Kawashima, K., Abdulla Daanaa, H. S., Yang, Z., & Akashi, H. (2025). Dinucleotide preferences underlie apparent codon preference reversals in the Drosophila melanogaster lineage. Proceedings of the National Academy of Sciences, 122(21), e2419696122. https://doi.org/10.1073/pnas.2419696122

Yamashita, H., Matsumoto, T., Kawashima, K., Daanaa, H. S. A., Yang, Z., & Akashi, H. (2024). Recent codon preference reversals in the Drosophila melanogaster lineage (p. 2024.10.10.617326). bioRxiv. https://doi.org/10.1101/2024.10.10.617326

Matsumoto, T., & Akashi, H. (2018). Distinguishing Among Evolutionary Forces Acting on Genome-Wide Base Composition: Computer Simulation Analysis of Approximate Methods for Inferring Site Frequency Spectra of Derived Mutations in Recombining Regions. G3, 8(5), 1755–1769. https://doi.org/10.1534/g3.117.300512

<!-- Links: below are not visible -->
[1]: https://doi.org/10.5281/zenodo.15274324

