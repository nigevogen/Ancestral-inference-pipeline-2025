# PPDA (Pipeline for Polymorphism and Divergence Analysis)

<!-- ![landing](tree.png) -->

This pipeline estimates ancestor states and counts genetic changes for individual type of mutation. __This package has a nice feature that estimates genetic changes for a single site on a particular branch so that we can relate genetic changes to external functional information such as protein binding affinity, ribosome stalling rate, nucleosome occupancy or so.__ A pipeline was initially developed to study evolution of codon usage, protein sequence and intron sequences in _Drosophila_ but it applicable to genomes for any kind of organisms.

## Processes

### I. Alignment preparation

An input file of aligned sequences is cleaned (remove gaps, Ns and other marked sites). Then segregating alleles in a population are collapsed to two variants to reduce information on allele frequency (allele frequency information is used later in BTW).

Example

__[INPUT] CDS alignment files__: FASTA formatted file of aligned coding sequences (CDS). Marker sequences can be included with suffix `_marker` in its sequence name.

```[Python]
>Dmel_sample1
ATGTCTGGGCGCGAGGGCGGTAAGAAGAAGCCTCTG
>Dmel_sample2
ATGTCCGGACGCGAGGGCGGCAAGAAGAAGCCTCTG
>Dmel_sample3
ATGTCTGGCCGCGAGGGCGGAAAGAAGAAGCCTCTG
>Dsim_sample1
ATGTCTGGACGCGAGGGCGGGAAGAAGAAGCCTCTG
>Dsim_sample2
ATGTCTGGACGCGAGGGCGGAAAGAAGAAGCCTCTG
>Dsim_sample3
ATGTCTGGACGCGAGGGCGGGAAGAAGAAGCCTCTG
>Dsim_sample4
ATGTCTGGACGCGAGGGCGGAAAG---AAGCCTCTG
>Dyak
ATGTCTGGACGCGAGGGCGGTAAAAAGAAGCCTCTG
>Dere
ATGTCTGGACGCNNNGGCGGTAAGAAGAAGCCTCTG
>conspos_marker
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCNNNNNN
```

Notes

- All stop codons have to be removed.
- A marker sequence has to include only `0` and `1`, or `N` and `C`. `N` is equivalent to `0` and sites aligned to any `0` or `N` are removed before ancestor inference.
- Names of sample sequences for a species have to start with the same prefix across alignment files.

After removing sites not used in ancestor inference (AI) and collapsing segregating alleles, the example alignment becomes like,

```[Python]
>Dyak
ATGTCTCGCGGCGGTAAAAAG
>Dere
ATGTCTCGCGGCGGTAAGAAG
>Dmel_collapse_a
ATGTCTCGCGGCGGCAAGAAG
>Dmel_collapse_b
ATGTCCCGCGGCGGTAAGAAG
>Dsim_collapse_a
ATGTCTCGCGGCGGAAAGAAG
>Dsim_collapse_b
ATGTCTCGCGGCGGGAAGAAG
```

At the second codon, two `TCT` and `TCC` of Dmel are collapsed to two codons of `TCT` and `TCC` in the output sequence. We shall call the output sequence __collapse sequence__ in this pipeline. Beside the second codons, the third codon was removed in the anlaysis, since there are more than two types of segregating alleles, which cannot be represented in collapse sequence. Although the 7th codon has four variants, since there are only two alleles segregating in each of population, the codon is kept in the output. The intermediate status of site filtering is not ouput by default, but you can see if site filtering was done properly by setting `True` to `verbose_output` option and temporally files are output.

Then the 3rd position of each codon is retained. The following ailgnment is used for the next ancestor inference program.

```[Python]
>Dyak
GTCCTAG
>Dere
GTCCTGG
>Dmel_collapse_a
GTCCCGG
>Dmel_collapse_b
GCCCTGG
>Dsim_collapse_a
GTCCAGG
>Dsim_collapse_b
GTCCGGG
```

__[INPUT] Alignment list__: A list of alignments to use from a given directory. Ancestor probabilities are estimated under assumption of parameter homogeneity across site and genes. In addition, this list is used to specify the order of concatenation of alignments, which is important to keep track of genomic positions of sites for which ancestor is estimated. An alignment list has to be formatted like,

```[Python]
1                           # Number of datasets (this should be always 1, currently)
prefix: Dmel_               # {prefix}{item}{suffix} should make an alignment file name.
suffix: .mpspye.aln         # From the first item, a file name "Dmel_13.mpspye.aln" is reconstructed.
>input_alignment 16      # >{dataset_name}\t{number of following items}
13                          # items are listed as follows.
14
15
1020
1023
1025
1028
1030
1031
1034
1036
1037
1040
1044
1047
1048
```

### II. Ancestor inference using BTW (Bifurcating Tree with Weighting) method

#### BASEML

Ancestor states are estimated for collapsed sequences using BASEML program in PAML (version 4).

__[INPUT] Sample tree__: A newick formatted file for collapse sequence. Nodes tagged with the same value are assumed to evolve under the shared parameters. This option is available at least from PAML version 4.9i. In this examle,

```[Python]
  6 1
((1 #0, 2 #1) #4, (3 #2, 4 #2) #5, (5 #3, 6 #3) #6) #7;
```

sample `3` and `4`, which represent Dmel collapse sequences, are assumed to evolve under the same evolutionary parameters but different from other branches. Similarly, Dsim samples are assumed to evolve under the same parameters. We used this option to restrict parameter range of population samples from other lineages. Please see [PAML manual][4] for details.

To run BASEML individually, you can use `run_BASEML.sh` script and please specify a path to temprate directory for BASEML run.

#### Weighting ancestor probabilities based on expected Site Frequecncy Spectra (SFS)

The ancestor probabilities inferred by BASEML were weighted by Site Frequency Spectrum (SFS) to taking allele frequency into account. This BTW method is proposed in [Matsumoto and Akashi, 2018][5] and Dr. Matsumoto has maintained the program since published. This method allows more accurate estimation of ancestor probabilities for population samples. Please see "BTW_instruction_TM_ver1.docx" for details.

To run BTW individually, please call `run_BTW_pipeline.sh` in terminal.

### III. Compute joint probabilities and Make Site Frequency Spectra

Mutations is counted by computing joint proability of across all scenarios for a site. You need to specify some parameters such as a path to sample sequence and information on population samples via a control file.

## Dependencies

- baseml (in [PAML][1] package)
- codonpaths
- pandas
- numpy

## Authors

- [Tomotaka Matsumoto][2]
- [Haruka Yamashita][3]

## References

- PAML: Yang Z., 2007 PAML 4: phylogenetic analysis by maximum likelihood. Mol. Biol. Evol. 24: 1586–1591. <https://doi.org/10.1093/molbev/msm088>
- BTW: Matsumoto T., Akashi H., 2018 Distinguishing Among Evolutionary Forces Acting on Genome-Wide Base Composition: Computer Simulation Analysis of Approximate Methods for Inferring Site Frequency Spectra of Derived Mutations in Recombining Regions. G3 (Bethesda): g3.300512.2017. <https://doi.org/10.1534/g3.117.300512>

[1]: http://abacus.gene.ucl.ac.uk/software/paml.html
[2]: https://github.com/tomotakamatsumoto
[3]: https://github.com/yamasampo
[4]: http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf
[5]: https://doi.org/10.1534/g3.117.300512
