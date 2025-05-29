# Ancestral-inference-pipeline-2025

A collection of codes used to infer lineage-specific nucleotide/codon changes from sequence alignments. The ancestral inference pipeline is designed for within- and between-species data. We provide one Python package and two pipelines that we used for analyses presented in Yamashita et al. 2025. The sequence alignment data used in the study is deposited at [Zenodo][1]. 

## Package / pipeline

This repository consists of one Python package, `refjoin`, and analysis pipeline, `ppda`. 

### refjoin

Combines two sequence alignments using a specified shared sequence. This package was developed by Kent Kawashima. `numpy` is required to use this package. We used this package to add within-species data to the reference sequence alignments. 

### ppda

Infers nucleotide changes on each branch of an assumed gene tree. This pipeline requires `BASEML` ([PAML website][2] and [Github][3]; Yang 2007). Ziheng Yang developed a new option that allows user-defined branches to share transition parameters. This new option is implemented in `BASEML` version 4.9 or later. BTW codes were developed by Tomotaka Matsumoto. See Matsumoto and Akashi 2018. This pipeline requires [alignmentrs][4] which allows fast processing of sequence alignments. `alignmentrs` was developed by Kent Kawashima. Other Python functions (e.g., for filtering aligned sites) and shell scripts to automate processes were developed by Haruka Yamashita. Pandas and numpy are also required to run this pipeline.

## References

Yamashita, H., Matsumoto, T., Kawashima, K., Abdulla Daanaa, H. S., Yang, Z., & Akashi, H. (2025). Dinucleotide preferences underlie apparent codon preference reversals in the Drosophila melanogaster lineage. Proceedings of the National Academy of Sciences, 122(21), e2419696122. https://doi.org/10.1073/pnas.2419696122

Yang, Z. (2007). PAML 4: Phylogenetic analysis by maximum likelihood. Molecular Biology and Evolution, 24(8), 1586–1591. https://doi.org/10.1093/molbev/msm088

Matsumoto, T., & Akashi, H. (2018). Distinguishing Among Evolutionary Forces Acting on Genome-Wide Base Composition: Computer Simulation Analysis of Approximate Methods for Inferring Site Frequency Spectra of Derived Mutations in Recombining Regions. G3, 8(5), 1755–1769. https://doi.org/10.1534/g3.117.300512

<!-- Links: below are not visible -->
[1]: https://doi.org/10.5281/zenodo.15274324
[2]: http://abacus.gene.ucl.ac.uk/software/paml.html
[3]: https://github.com/abacus-gene/paml
[4]: https://github.com/kentwait/alignmentrs

