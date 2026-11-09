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
Please see the project home page for usage details and installation instructions: https://github.com/Ecogenomics/CheckM/wiki

We do not recommend installing CheckM from the master branch. This may be unstable. Please install an official release of CheckM or use pip.

### Dependencies

CheckM requires [HMMER](http://hmmer.org/) (>=3.1), [Prodigal](https://github.com/hyattpd/Prodigal) (2.6.3) and, from CheckM 1.3.0, [gPlace](https://github.com/cmc-aau/gPlace) (>=0.1.0) to be on your system path. gPlace replaces [pplacer](https://github.com/matsen/pplacer), which was required by earlier versions of CheckM: it gives the same placements as pplacer, bit for bit, using less memory and time. gPlace requires Python >=3.12 and can be installed with `pip install gtdb-gplace`; it can be installed in a separate environment from CheckM as long as the `gplace` command is on your system path.

CheckM >= 1.3.0 requires Python >= 3.11 and [xopen](https://github.com/pycompression/xopen) (>= 2.0), which is installed with CheckM. From CheckM 1.3.0, most output files are gzip-compressed (e.g., `lineage.ms.gz`, `bins/<genome_id>/genes.faa.gz` and `storage/bin_stats_ext.tsv.gz`); CheckM commands still accept the uncompressed output of earlier versions and the uncompressed file names (e.g., `checkm qa ./output/lineage.ms ./output`).

## Estimating quality of CPR genomes

Information about obtaining improved quality estimates for CPR (Patescibacteria) genomes can be found here:
https://github.com/Ecogenomics/CheckM/wiki/Workflows#using-cpr-marker-set

## Migration to Python 3

CheckM has been ported to Python 3 to accomodate Python 2 reaching [end of life](https://pythonclock.org/) on January 1, 2020. CheckM >=1.1.0 requires Python 3. Python 2 will no longer be actively supported. Apologies for any issues this may cause.

Massive thanks to [baudrly](https://github.com/baudrly), [Vini Salazar](https://github.com/vinisalazar), and [Asaf Peer](https://github.com/asafpr) for initial Python 2 to 3 porting.

### Python 2 to 3 Validation

Porting of CheckM to Python 3 was validation on a set of 1,000 genomes randomly select from the [GTDB](https://gtdb.ecogenomic.org/) R89 representative genomes. Results were compared to those generated with CheckM v1.0.18, the last Python 2 version of CheckM. Identical results were obtained for the 'lineage_wf', 'taxonomy_wf', and 'ssu_finder' methods across this set of test genomes. Other CheckM methods have been executed on a small set of 3 genomes to verify they run to completion under Python 3. 

### Removed Functionality

The following features have been removed from CheckM v1.1.x in order to simplify the code base and focus CheckM and support requests on critical functionality:
 * bin_qa_plot: non-critical, rarely used plot which does not scale to the large numbers of MAGs now being recovered
 * par_plot: non-critical plot and the same information is better presented in the reference distribution plots
 * cov_pca, tetra_pca: alternatives to these static plots exist in tools such as [Anvi'o](http://merenlab.org/software/anvio/)
 * len_plot: rarely used plot which is largely redundant with the len_hist and nx_plot plots
  * bin_union, bin_compare: feature rich alternative now exist such as [DAS Tool](https://github.com/cmks/DAS_Tool) and [UniteM](https://github.com/dparks1134/UniteM)

## Bug Reports

Please report bugs through the GitHub issues system. 

Copyright © 2014 Donovan Parks, Connor Skennerton, Michael Imelfort. See LICENSE for further details.
