# CheckM

[![version status](https://img.shields.io/pypi/v/checkm-genome.svg)](https://pypi.python.org/pypi/checkm-genome)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/checkm-genome.svg?color=green)](https://anaconda.org/bioconda/checkm-genome)
[![Downloads](https://pepy.tech/badge/checkm-genome/month)](https://pepy.tech/project/checkm-genome)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/checkm-genome.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/checkm-genome)

## Installing and using CheckM
Please see the project home page for usage details and installation instructions: https://github.com/Ecogenomics/CheckM/wiki

We do not recommend installing CheckM from the master branch. This may be unstable. Please install an official release of CheckM or use pip.

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
