# CheckM for Python 3

This branch is the initial release of CheckM for Python 3. CheckM has been ported to Python 3 as Python 2 is reaching [end of life](https://pythonclock.org/) on January 1, 2020. CheckM will be updated from v1.0.18 to v1.1.x on Dec. 1, 2019. CheckM v1.1.x will require Python 3 and there are no plans to support Python 2 moving forward. Apologies for any issues this may cause.

Massive thanks to [baudrly](https://github.com/baudrly), [Vini Salazar](https://github.com/vinisalazar), and [Asaf Peer] (https://github.com/asafpr) for initial Python 2 to 3 porting.

## Python 2 to 3 Validation

Porting of CheckM to Python 3 was validation on a set of 1,000 genomes randomly select from the [GTDB](https://gtdb.ecogenomic.org/) R89 representative genomes. Results were compared to those generated with CheckM v1.0.18, the last Python 2 version of CheckM. Identical results were obtained for the 'lineage_wf', 'taxonomy_wf', and 'ssu_finder' methods across this set of test genomes. Other CheckM methods have been executed on a small set of 3 genomes to verify they run to completion under Python 3. 

## Removed Functionality

The following features will be removed from CheckM v1.1.x in order to simplify the code base and help focus the application and support requests on critical functionality:
 * bin_qa_plot: non-critical, rarely used plot which does not scale to the large numbers of MAGs now being recovered
 * par_plot: non-critical plot and the same information is better presented in the reference distribution plots
 * cov_pca, tetra_pca: alternatives to these static plots exist in tools such as [Anvi'o](http://merenlab.org/software/anvio/)
 * len_plot: rarely used plot which is largely redundant with the len_hist and nx_plot plots
  * bin_union, bin_compare: feature rich alternative now exist such as [DAS Tool](https://github.com/cmks/DAS_Tool) and [UniteM](https://github.com/dparks1134/UniteM)

## Bug Reports

Please report any bugs or observed discrpencies between this branch and CheckM v1.0.18 as a GitHub issue. 

Copyright © 2014 Donovan Parks, Connor Skennerton, Michael Imelfort. See LICENSE for further details.
