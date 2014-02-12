#Overview

Determining if contigs 'binned' into a putative genome is incomplete or contains sequence data from other genomes is an open problem in metagenomic assembly. These issues must be addressed to make informed inferences about the gene content and metabolic function of putative population genomes. CheckM provides a set of tools for assessing the quality of putative genome bins. It provides robust estimates of genome completeness and contamination using lineage-specific sets of single-copy, ubiquitious genes and explicitly accounting for gene co-location. Assessment of bin quality can also be examined using plots depicting key genomic characteristics (e.g., GC, coding density) which highlight contigs outside the expected distributions of a typical genome. CheckM also provides tools for identifying bins that are likely candidates for merging based on marker set compatibility, similarity in genomic characteristics, and proximity within a reference genome tree.

#Installation and Dependancies

#Running Instructions

# Marker Sets

## Generation of taxon specific marker sets
The taxon specific marker sets that are used by CheckM were generated
using the metadata from the IMG 4.0 release.  First all of the complete
genomes were separated based on the taxonomy string found in the
metadata.  Next the prevalence of pfam domains for each taxon with at
least 5 genomes was calculated todetermine which were found in at least
95% of all genomes.  This list of domainswas further segregated based on
whether that pfam was strictly single copy in each genome or not.

### Caveats

## Using Taxon Specific Marker Sets

# Workflow


#Licence and referencing

Project home page, info on the source tree, documentation, issues and
how to contribute, see https://github.com/dparks1134/CheckM

This software is currently unpublished. Please contact me at
donovan_parks_at_gmail_dot_com for more information about
referencing this software.

Copyright © 2013 Donovan Parks, Connor Skennerton, Michael Imelfort. See LICENSE.txt
for further details.
