# Overview
Determining if contigs 'binned' into a putative genome is incomplete or contains sequence data from other genomes is an open problem in metagenomic assembly. These issues must be addressed to make informed inferences about the gene content and metabolic function of putative population genomes. CheckM provides a set of tools for assessing the quality of putative genome bins. It provides robust estimates of genome completeness and contamination using lineage-specific sets of single-copy, ubiquitious genes and explicitly accounting for gene co-location. Assessment of bin quality can also be examined using plots depicting key genomic characteristics (e.g., GC, coding density) which highlight contigs outside the expected distributions of a typical genome. CheckM also provides tools for identifying bins that are likely candidates for merging based on marker set compatibility, similarity in genomic characteristics, and proximity within a reference genome tree.

# Installation and Dependancies
CheckM has been developed to work on a Linux system. It may work elsewhere, but I haven't tried. If you try it somewhere else then please let me know. I'd like to keep this list up-to-date.

### Using Pip

This is the recommended method for obtaining CheckM as it will automaticaly install a number of dependencies. You'll need to have pip installed on your system. After that, open a terminal and it is as easy as typing:

$ sudo pip install checkm-ACE

CheckM makes use of the following 3rd party programs which must be available on the system path:

    * hmmer >= 3.0
    * prodigal >= 2.6
    * pplacer >= 1.1

### Installing from source

If you prefer this type of thing you can always try install from source directly. You will need the following dependencies:

    * numpy >= 1.8.0
    * scipy >= 0.9.0
    * matplotlib >= 1.1.0
    * dendropy >= 3.12.0
    * pysam >= 0.7.4

The 3rd party dependencies listed above must also be install and placed on the system path. Clone the repo from github:

$ git clone https://github.com/dparks1134/CheckM.git

Then change into the CheckM directory and type:

$ sudo python setup.py install

# Using CheckM

You can checkout our new and expanding manual (currently under construction), or read on if you're in a rush.
Before using CheckM you need to generate putative genome bins from your metagenomic dataset(s). Our companion tool GroopM can be used for this task.CheckM assumesgenome binsare contained in a single folder with a common file extension.

### Standard workflow

CheckM was developed to be highly flexible in the marker set used for assessing bin completeness and contamination. Functionality is provided for calculatinglineage- and taxonomic-specific sets of marker, or you can supply your own HMM file specifying a custom set of markers.The standard workflow for using lineage-specific markers is as follows:

    * tree - infer the position of each bin within a reference genome tree
    * tree_qa - examine the placementof bins within the genome tree and their taxonomic affliations
    * lineage_set - infer a lineage-specific marker set for each bin
    * analyze - identify markers in each bin
    * qa - assess bin completeness and contamination

For more information on these commands type:

$ checkm OPTION -h

### Plots

A suite of plots is provided for investigating the quality of genome bins. We recommend that at least the following plots be considered:

    * bin_qa_plot - illustratescompleteness, contamination, and strain heterogeneity of each bin
    * dist_plot - produces a plot for each bin illustrating the bin's GC, coding density, and tetranucleotide signature distributions relative to expected distributions as determined from a set of reference genomes

### Bin mergers

Automated methods for binning contigs tend to be conservative which can result in a genome being split across multiple bins. To help identify bins that are likely candidates formerging, CheckM provides functionality for identifying bins with complementary sets of markers:

    * merge - identify bins with complementary markers

The merge command requiresbins to be assessed using a common set of marker genes. Typically, this will be done with with both a BacterialandArchaeal set generated with the taxon_set command (see manual for further details).

### Going further
CheckM provides a number of additional plots and functionality for performing tasks such as bin refine, identifying unbinned contigs, and producing community profiles. See the manual for more details.

# Licence and referencing

Project home page, info on the source tree, documentation, issues and
how to contribute, see https://github.com/dparks1134/CheckM

This software is currently unpublished. Please contact me at
*donovan_parks_at_gmail_dot_com* for more information about
referencing this software.

Copyright © 2014 Donovan Parks, Connor Skennerton, Michael Imelfort. See LICENSE.txt
for further details.
