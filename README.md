[IN ACTIVE DEVELOPMENT AND NOT CURRENTLY RECOMMEND FOR USE!]

# CheckM for Python 3

This branch is the initial release of CheckM for Python 3. CheckM has been ported to Python 3 as Python 2 is reaching [end of life](https://pythonclock.org/) on January 1, 2020. CheckM will be updated from v1.0.18 to v1.1.0 on Dec. 1, 2019. CheckM v1.1.0 will require Python 3 and there are no plans to maintain a Python 2 compatible version of CheckM. Apologies for any issues this may cause.

Massive thanks to baudrly, vinisalazar, and asafpr for initial Python 2 to 3 porting.


## Python 2 to 3 Validation

Porting of CheckM to Python 3 was validation on a set of 1,000 genomes randomly select from the [GTDB](https://gtdb.ecogenomic.org/) R89 representative genomes. Results were compared to those generated with CheckM v1.0.18, the last Python 2 version of CheckM. Identical results were obtained for the 'lineage_wf', 'taxonomy_wf', and 'ssu_finder' methods across this set of test genomes. Other CheckM methods have been executed on a small set of 3 genomes to verify they run to completion under Python 3. 

## Bug Reports

Please report any bugs or observed discrpencies between this branch and CheckM v1.0.18 as a GitHub issue. 

Copyright © 2014 Donovan Parks, Connor Skennerton, Michael Imelfort. See LICENSE for further details.
