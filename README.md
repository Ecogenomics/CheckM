#Overview

One of the problems faced in metagenomic assembly is the uncertainty
that collections of contigs that have been 'binned' into putative
genomes may be contaminated from another organism or represent an
incomplete genome.  Both are major problems that need to be addressed
before inferences about metabolic function can be made by the presence
or absence of genes. CheckM, is a set of tools that can assess the
contamination and completeness of genome bins using conserved genes found
in a single copy for different taxonomic lineages.  More closely
related bacteria and archaea share a higher proportion of their gene
complement and therefore estimations of completeness based on these
markers will be more accurate than estimates based on gene markers
that are universal.

#Installation and Dependancies

CheckM requires both `biopython` and `simpleHMMER` python packages be
intalled as well as `hmmer3` binaries to be in your `PATH`.  With those
dependancies filling download the git repository to your system and run
```
python setup.py install
```

#Running Instructions

CheckM has a single executable `checkM` that has a number of
subcommands.  There are three subcommands that the use will use
regularly:
* `build` - Run `hmmer` of a marker set against genome(s) of interest
* `qa` - Filter `hmmer` results and output a summary
* `align` - Align identified markers against the HMM model
CheckM can be run in two different ways.  First is for the user to
supply a `hmmer3` compatible model file of markers.  The second way is
to build a database of markers for all taxonomic lineages and use those
markers for accurate estimates of completeness and contamination.  By
default CheckM prefers to use the later of the two options.


## Setting up the Database
`checkM` has a subcommand called `makeDB` that will take freely
available genomic data and calculate the conserved markers for each
taxonomic lineage.  To generate the database you will need to download
the genomes from IMG and have the Pfam-A domains as well.

Assuming that you have either made a database you should set the path
with the environmental variable `CHECKM_DB`.  To check whether this is
already set, try the following command:
```
$ echo $CHECKM_DB
```
If you see a file path then the CheckM database is already set up on
your system.  Otherwise you can setup the database using the following:
```
export CHECKM_DB=/full/path/to/your/checkm_database
```
To make that setting perminant add the above line into your `.bashrc` file.

# Marker Sets
Marker sets are a list of genes/protein domains that are found in all
genomes for a particular taxonomic rank.  Some marker sets may be
universal or be specific to a family or genus.  There are two marker
sets that get used routinely; one from
[phylosift](http://phylosift.wordpress.com) that contains 37 markers
that are universal, and a list compiled by [Dupont et al.
(2012)](http://www.ncbi.nlm.nih.gov/pubmed/?term=22170421) which
contains 111 markers found in the majority of Bacteria. A more accurate
view of contamination and completeness of a genome bin may be achieved
by using markers that are targeted to a particular microbial lineage.


## Generation of taxon specific marker sets
The taxon specific marker sets that are used by CheckM were generated
using the metadata from the IMG 4.0 release.  First all of the complete
genomes were separated based on the taxonomy string found in the
metadata.  Next the prevalence of pfam domains for each taxon with at
least 5 genomes was calculated todetermine which were found in at least
95% of all genomes.  This list of domainswas further segregated based on
whether that pfam was strictly single copy in each genome or not.

### Caveats
  * These lists use pfams, which means that a marker may not represent a whole
    gene.  For example there may be two markers, one from the C-terminal and
    the other from the N-terminal domain of the same protein.
  * using a minimum of 5 genomes for a taxon is probably too low.  Be careful
    about how specific a marker set you use.  Think about the total number of
    genomes that the list was constructed against, if in doubt use a Domain
    or phylum level marker set.
  * based solely on the genomes found in IMG 4.0
  * I haven't done any annotations myself, I'm relying solely on the annotation
    pipeline in IMG
  * The taxonomy is also straight from IMG and this may be inaccurate (maybe?),
    which would throw off all of the calculations

## Using Taxon Specific Marker Sets
Taxon specific marker sets can be invoked by using the `-T` option on
the command-line.  The taxonomy must be the full taxonomy in Greengenes
format, for example: `k__Bacteria;p__Proteobacteria`.  Spaces are also
allowed between taxonomic ranks: `k__Bacteria; p__Proteobacteria`. To
get a universal marker set (common to all bacteria and archaea) use the
special taxonomy 'universal' as the argument to `-T`
NOTE: Do **not** specify a taxonomy that is lower than the genus level!

# Workflow
This is a workflow that I use for determining the completeness and
contamination of all of my genome bins.  First start by using one of the
universal marker sets such as the markers used by phylosift or by
specifying the taxonomy as 'universal' in the `-T` flag of `checkM.
```
checkM build -T universal <OUTPUT_FOLDER> genome.fna

checkM qa -T universal <OUTPUT>
```

From the output of the `qa` command you should be able to see which bins
are 'near complete'.  As a general rule I only look genome bins that are
at least 80% complete and less than 10% contamination.  From this point
I make a tree using the concatenation of the  markers and determine the
putative taxonomy of my genome bins.

For example lets assume that one of your genome bins is a
__Rhodocyclacae__ based on the taxonomy in the tree.  Now we can use the
specific single copy genes for __Rhodocyclacae__ to get more accurate
statistics on completeness and contamination. (There are 428
__Rhodocyclacae__ specific single copy genes, verses the 39 universal)
```
checkM build -T 'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Rhodocyclales;f__Rhodocyclaceae' <OUTPUT_FOLDER> genome.fna

checkM qa -T 'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Rhodocyclales;f__Rhodocyclaceae' <OUTPUT_FOLDER>
```
This output should be a more accurate view of the contamination and completeness of your genome.

# Manual
## Build
```
checkM build [-h] [-T TAXONOMY] [-d DATABASE] [-H HMM] [-b BIN_FOLDER]
                    [-x EXTENSION] [-t THREADS] [-p PREFIX] [-v]
                    out_folder [FASTA [FASTA ...]]

Parse binned contigs

positional arguments:
  out_folder            a place to write output files
  FASTA                 A list of fasta formatted files to be used as the
                        input (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -T TAXONOMY, --tax-string TAXONOMY
                        Specify a Greengenes taxonomy string for the lineage
                        of interest (default: None)
  -d DATABASE, --database DATABASE
                        Specify a different path for the marker database. The
                        default is to use the value of the environmental
                        variable "CHECKM_DB" (default: None)
  -H HMM, --hmm HMM     hmm used in the search (default: None)
  -b BIN_FOLDER, --bin_folder BIN_FOLDER
                        folder containing bins to check (fasta formatted)
                        (default: None)
  -x EXTENSION, --extension EXTENSION
                        used in conjuntion with -b to specify the extension
                        used on files to be parsed (default: fa)
  -t THREADS, --threads THREADS
                        max number of active threads (default: 1)
  -p PREFIX, --prefix PREFIX
                        prefix used for naming output files (default: )
  -v, --verbose         print more (default: False)
```

## QA
```
checkM qa [-h] [-o OUT_FORMAT] [-T TAXONOMY] [-d DATABASE] [-H HMM]
                 [-p PREFIX] [-e E_VALUE] [-l LENGTH] [-f FILE] [-v]
                 out_folder

Do QA on pre-processed contigs

positional arguments:
  out_folder            folder specified during build command

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_FORMAT, --output-format OUT_FORMAT
                        Change the output format to one of the following:
                        1. Simple Summary showing bin name, counts, completeness
                           and contamination
                        2. List of marker names and their counts
                        3. Matrix of marker counts
                        4. Tabular list of bin name, marker name, contig name
                        5. Tabular list showing contigs that contain more than
                           one copy of the same marker. Format: bin name,
                           contig name, {Marker name, count}... (default: 1)
  -T TAXONOMY, --tax-string TAXONOMY
                        Specify a Greengenes taxonomy string for the lineage
                        of interest (default: None)
  -d DATABASE, --database DATABASE
                        Specify a different path for the marker database. The
                        default is to use the value of the environmental
                        variable "CHECKM_DB" (default: None)
  -H HMM, --hmm HMM     hmm used in the search (default: None)
  -p PREFIX, --prefix PREFIX
                        prefix which was used for naming output files
                        (default: )
  -e E_VALUE, --e_value E_VALUE
                        e value cut off (default: 1e-10)
  -l LENGTH, --length LENGTH
                        percent overlap between target and query (default:
                        0.7)
  -f FILE, --file FILE  print results to file (default: STDOUT)
  -v, --verbose         print more (default: False)
```
### Output Format Examples
#### One
#### Two
#### Three
#### Four
#### Five

## Align
```
checkM align [-h] [-T TAXONOMY] [-d DATABASE] [-H HMM] [-p PREFIX]
                    [-e E_VALUE] [-l LENGTH] [-v] [-s] [-o OUT_FORMAT] [-c]
                    [-b]
                    out_folder

Create alignments of hmms to the matched markers

positional arguments:
  out_folder            folder specified during build command

optional arguments:
  -h, --help            show this help message and exit
  -T TAXONOMY, --tax-string TAXONOMY
                        Specify a Greengenes taxonomy string for the lineage
                        of interest (default: None)
  -d DATABASE, --database DATABASE
                        Specify a different path for the marker database. The
                        default is to use the value of the environmental
                        variable "CHECKM_DB" (default: None)
  -H HMM, --hmm HMM     hmm used in the search (default: None)
  -p PREFIX, --prefix PREFIX
                        prefix which was used for naming output files
                        (default: )
  -e E_VALUE, --e_value E_VALUE
                        e value cut off (default: 1e-10)
  -l LENGTH, --length LENGTH
                        percent overlap between target and query (default:
                        0.7)
  -v, --verbose         print more (default: False)
  -s, --separate-files  Output a separate file for each of the hits to a HMM
                        model. (default: False)
  -o OUT_FORMAT, --output-format OUT_FORMAT
                        The output format passed to hmmalign (default:
                        PSIBLAST)
  -c, --include-consensus
                        Include the HMM consensus sequence (generated by
                        hmmemit) in the alignment (default: False)
  -b, --best-hit        Only align the best match to the HMM even if there are
                        multiple matches that pass the threshold (default:
                        False)
```

## makeDB
```
checkM makeDB [-h] [-f] [-m] [-i IMG_BASE] [-M METADATA]
                     [--minimum-genome-count MIN_GENOME]
                     [--minimum-taxon-conservation MIN_CONS] [-p PFAM]
                     [-o DATABASE]

optional arguments:
  -h, --help            show this help message and exit
  -f, --force           overwite markers or taxons ifthey are already stored
                        in the database (default: False)
  -m, --make_tables     make new tables from scratch (default: False)
  -i IMG_BASE, --input-prefix IMG_BASE
                        The root directory to all the IMG genomes (default:
                        None)
  -M METADATA, --metadata METADATA
                        The IMG metadata file (default: None)
  --minimum-genome-count MIN_GENOME
                        A taxonomic lineage must have at least this number of
                        genomes for a marker set to be made (default: 5)
  --minimum-taxon-conservation MIN_CONS
                        Percent of genomes in taxon that must contain the
                        marker (default: 0.95)
  -p PFAM, --pfam-database PFAM
                        path to Pfam-A file (default: None)
  -o DATABASE, --output-database DATABASE
                        output sqlite3 database name (default: None)
```


#Licence and referencing

Project home page, info on the source tree, documentation, issues and
how to contribute, see http://github.com/Ecogenomics/CheckM

This software is currently unpublished. Please contact me at
donovan_parks_at_gmail_dot_com for more information about
referencing this software.

Copyright © 2012, 2013 Donovan Parks, Michael Imelfort, Connor Skennerton. See LICENSE.txt
for further details.
