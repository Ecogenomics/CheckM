# CheckM

[![Unmaintained](http://unmaintained.tech/badge.svg)](http://unmaintained.tech/)   
[![version status](https://img.shields.io/pypi/v/checkm-genome.svg?cacheSeconds=3600)](https://pypi.python.org/pypi/checkm-genome)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/checkm-genome.svg?color=green&cacheSeconds=3600)](https://anaconda.org/bioconda/checkm-genome)
[![Downloads](https://pepy.tech/badge/checkm-genome/month)](https://pepy.tech/project/checkm-genome)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/checkm-genome.svg?style=flag&label=BioConda%20install&cacheSeconds=3600)](https://anaconda.org/bioconda/checkm-genome)

## CheckM v2 has been released

CheckM v2 has been released and can be found [here](https://github.com/chklovski/CheckM2). CheckM v1 (this repository) estimates genome completeness and contamination based on the presence or absence of marker genes, i.e., genes that are typically ubiquitous and single copy. In contrast, [CheckM v2](https://github.com/chklovski/CheckM2) uses machine learning models to estimate completeness and contamination. Our [assessment](https://www.nature.com/articles/s41592-023-01940-w) indicates CheckM v2 will be more accurate in general. However, you may consider running both CheckM v1 and v2 in some cases since these are independent methodologies for estimating genome quality. It can also be insightful to explore the marker genes that are missing in your genomes, especially if you have multiple related genomes and are interested in systematic loss of an otherwise ubiquitious gene. 

Outside of minor updates to fix compatibility issues with external dependencies no maintanence is intended for CheckM v1.

## Installing and using CheckM
Please see the CheckM documentation for usage details and installation instructions: https://ecogenomics.github.io/CheckM

We do not recommend installing CheckM from the master branch. This may be unstable. Please install an official release of CheckM or use pip.

### Dependencies

CheckM requires [HMMER](http://hmmer.org/) (>=3.1), [Prodigal](https://github.com/hyattpd/Prodigal) (2.6.3) and, from CheckM 1.3.0, [gPlace](https://github.com/cmc-aau/gPlace) (>=0.1.0) to be on your system path. gPlace replaces [pplacer](https://github.com/matsen/pplacer), which was required by earlier versions of CheckM: it gives the same placements as pplacer, bit for bit, using less memory and time. gPlace requires Python >=3.12 and can be installed with `pip install gtdb-gplace`; it can be installed in a separate environment from CheckM as long as the `gplace` command is on your system path.

CheckM >= 1.3.0 requires Python >= 3.11 and [xopen](https://github.com/pycompression/xopen) (>= 2.0), which is installed with CheckM. From CheckM 1.3.0, most output files are gzip-compressed (e.g., `lineage.ms.gz`, `bins/<genome_id>/genes.faa.gz` and `storage/bin_stats_ext.tsv.gz`); CheckM commands still accept the uncompressed output of earlier versions and the uncompressed file names (e.g., `checkm qa ./output/lineage.ms ./output`).

## Estimating quality of CPR genomes

Information about obtaining improved quality estimates for CPR (Patescibacteria) genomes can be found here:
https://ecogenomics.github.io/CheckM/workflows.html#using-the-cpr-marker-set

## Bug Reports

Please report bugs through the GitHub issues system. 

Copyright © 2014 Donovan Parks, Connor Skennerton, Michael Imelfort. See LICENSE for further details.
