# Change log

All published releases of CheckM are listed below, most recent first. Dates are those of the [GitHub release](https://github.com/Ecogenomics/CheckM/releases) or, where there is none, of the upload to [PyPI](https://pypi.org/project/checkm-genome/#history). Some version numbers were not released separately; their changes are listed under the release that first included them.

## 1.4.0 (2026-11-23)

HMMER is now run through [pyhmmer](https://github.com/althonos/pyhmmer), so the HMMER programs no longer need to be installed. Results are unchanged: `quality_report.tsv` of `lineage_wf` is identical to CheckM 1.3.2 on 128 and 1,000 test genomes, and that of `taxonomy_wf` is identical on 128 genomes. On 1,000 genomes with 64 threads, `lineage_wf` takes 9% less time than CheckM 1.3.2 (17.0 instead of 18.7 min). See [Migration to pyhmmer](https://ecogenomics.github.io/CheckM/pyhmmer_migration.html) and the [performance results](https://ecogenomics.github.io/CheckM/performance.html).

- Adds pyhmmer (>= 0.12.0) as a dependency, which is installed with CheckM and includes HMMER 3.4. It is available for Linux and macOS. `hmmsearch`, `hmmalign`, `hmmfetch` and `nhmmer` are no longer required on the system path.
- Identifying marker genes takes less time but needs more memory. With 64 threads and 1,000 genomes, this step of `analyze` takes 11% less time (11.6 instead of 13.1 min) and needs 5.9 GiB instead of 3.1 GiB. The peak memory of `tree` and `lineage_wf` is still set by gPlace, unless its memory is limited with `--gplace_max_memory` (e.g., 6.3 GiB instead of 5.3 GiB for `lineage_wf` with `--gplace_max_memory 6G`).
- gPlace >= 0.1.2 is required (previously >= 0.1.0). Its memory no longer differs between runs on the same genomes (13.6 GiB with 64 threads and 1,000 genomes; gPlace 0.1.0 used 14.3 to 18.4 GiB) and it is faster when its memory is limited (19 s instead of about 1 min with `--gplace_max_memory 6G`).
- E-values that `hmmsearch` prints as `0` may be written as very small numbers (e.g., `5.2e-313`) in the HMMER tables (e.g., `bins/<genome_id>/hmmer.analyze.txt.gz`). All other values are unchanged.
- The report written with `--ali` (e.g., `bins/<genome_id>/hmmer.analyze.ali.txt.gz`) is now written by CheckM in the layout of the output of `hmmsearch`. It has the same hits and alignments, but not the header and search statistics of `hmmsearch`.
- `ssu_finder` runs nhmmer through pyhmmer and writes the hit tables itself. `ssu_summary.tsv` and `ssu.fna` are identical to CheckM 1.3.2 on 128 genomes. It is slightly slower (15 s instead of 10 s with 64 threads).
- The Docker image no longer installs HMMER.

## 1.3.2 (2026-11-19)

Reduces the memory needed to call genes and fixes occasional high memory usage by `analyze`. Results are unchanged: all outputs of `lineage_wf` and `taxonomy_wf` are identical to CheckM 1.3.1 on 128 test genomes, and `quality_report.tsv` is identical on 1,000 genomes.

- Gene calling (the first step of `tree`, `lineage_wf` and `taxonomy_wf`) no longer starts the largest genomes all at once. Genomes are still processed from largest to smallest, but a genome is only started if the estimated memory of Prodigal fits within a budget set from the sizes of the genomes, and smaller genomes are processed in the meantime. With 64 threads and 1,000 genomes, gene calling needs 5.5 GiB instead of 8.3 GiB and takes 4% longer (about 10 s).
- `analyze` (including the `analyze` step of `lineage_wf` and `taxonomy_wf`) could occasionally need much more memory: its worker processes sometimes made private copies of the memory of the main process when they ran Python's garbage collector, as fixed for other steps in 1.3.1 (e.g., 10.3 GiB instead of ~3 GiB in one run on 1,000 genomes with 64 threads). Worker processes now keep sharing this memory.
- Coding density is calculated from the positions of genes instead of an array with 8 bytes per base of each contig, which reduces the memory used by each thread when calling genes and calculating genome statistics (e.g., by ~0.07 GiB for a 10 Mbp genome).

## 1.3.1 (2026-11-19)

Fixes high memory usage with many genomes and a failure with many novel genomes. Results are unchanged: all outputs of `lineage_wf`, `taxonomy_wf`, every `qa` output format and option, and `tree_qa` are identical to CheckM 1.3.0 on 128 test genomes, and the `qa` results and multi-copy marker alignments are identical on 10,000 genomes.

- Less memory is needed by `analyze` (the "Aligning marker genes with multiple hits in a single bin" and genome statistics steps), `qa` and the `qa` step of `lineage_wf` and `taxonomy_wf`. Previously, the HMM hits of all genomes were read before the work was split across threads, and each thread made its own copy of this data, so memory grew with the number of genomes multiplied by the number of threads. HMM hits are now read and reported one genome at a time, marker sets and HMM models are shared between genomes, and threads no longer copy the data of the main process. With 64 threads and 10,000 genomes, the multi-hit alignment step needs 3.9 GiB instead of about 40 GiB and `qa` needs 3.9 GiB instead of 29.5 GiB. With 32 threads and 1,000 genomes, the multi-hit alignment step needs 3.3 GiB instead of 35 GiB.
- `lineage_set` and `lineage_wf` no longer fail with a `RecursionError` when many genomes are placed on the same branch of the reference tree (e.g., 10,000 novel genomes gave a tree nested ~1,000 levels deep).

## 1.3.0 (2026-11-10)

**CheckM 1.3.0 changes some arguments and outputs, so scripts that run CheckM or read its output files may need to be updated.** Completeness and contamination estimates are identical to CheckM 1.2.8 to 1.2.10 on all 1,128 test genomes. On 1,000 genomes with 64 threads, `lineage_wf` takes 19% less time and needs 61% less memory (47 to 18 GiB) and 72% less disk space (4.4 to 1.2 GB). See the [performance results](https://ecogenomics.github.io/CheckM/performance.html) and the [quality estimate changes](https://ecogenomics.github.io/CheckM/quality_estimate_changes.html) compared with CheckM 1.2.7 and earlier.

- pplacer has been replaced by [gPlace](https://github.com/cmc-aau/gPlace) (PyPI: `gtdb-gplace`), which gives the same placements as pplacer, bit for bit, using less memory and time. gPlace >= 0.1.0 must be on the system path (it requires Python >= 3.12 and can be installed in a separate environment), and pplacer and guppy are no longer required.
- The `--pplacer_threads` option of the `tree` and `lineage_wf` commands has been removed. gPlace uses the number of threads given by `-t`, without the previous limit of 64 threads.
- The names of gPlace's output files are those of pplacer, apart from compression (e.g., `storage/tree/concatenated.pplacer.json.gz`; see below).
- The `--reduced_tree` option of the `tree` and `lineage_wf` commands has been removed. Bins are always placed in the full reference genome tree, which gPlace does with about 18 GiB of memory with 64 threads, mostly for the reference tree (pplacer needed about 45 GiB with the full tree and less than 15 GiB with the reduced tree). Use `--gplace_max_memory` to place bins with less memory (e.g., `8G`), which gives identical results.
- New `--gplace_max_memory` option for the `tree` and `lineage_wf` commands sets the memory gPlace may use (e.g., `64G`). If the placement needs more, gPlace places bins tile by tile within this limit, which gives identical results.
- `checkm test` runs gPlace with `--gplace_max_memory 6G`, so it needs about 5 GiB of memory.
- When gPlace fails, its error message (e.g., that `--gplace_max_memory` is too small) is included in CheckM's error message.
- Progress is reported with tqdm progress bars (e.g., `Processing bins: 100%|##########| 128/128 [01:15<00:00, 1.70bin/s]`), which adds tqdm as a dependency. Progress is no longer reported when console output is suppressed with `--quiet`.
- Most output files are gzip-compressed (compression level 1, using xopen and ISA-L), which reduces the size of the output directory ~3.5x (e.g., from 596 MB to 169 MB for 128 genomes) and has little effect on run time. Results are unchanged. The compressed files are:
  - `bins/<genome_id>/*` (e.g., `genes.faa.gz`, `genes.gff.gz`, `hmmer.tree.txt.gz`, `hmmer.analyze.txt.gz`)
  - `lineage.ms` (now `lineage.ms.gz`) and the marker files of `lineage_set` and `taxon_set` (`.gz` is appended to the given name if not present)
  - `storage/*.tsv` (e.g., `storage/bin_stats_ext.tsv.gz`)
  - all files in `storage/tree` except the `*.log` and `*.out` files (e.g., `storage/tree/concatenated.tre.gz`)
  - the `storage/aai_qa` directory, which is replaced by the archive `storage/aai_qa.tar.gz`
- All commands still read the uncompressed output of earlier versions of CheckM and accept the uncompressed file names (e.g., `checkm qa ./output/lineage.ms ./output` reads `lineage.ms.gz`; if both `lineage.ms` and `lineage.ms.gz` exist, the most recently modified file is used). Scripts that read these files directly by their previous names must be updated.
- Adds xopen (>= 2.0) as a dependency and requires Python >= 3.11.
- `lineage_wf` and `taxonomy_wf` write their results to `<output_dir>/quality_report.tsv`, the tab-separated table produced by `qa -o 2 --tab_table`, and no longer print a results table to the console. Their `-f/--file` and `--tab_table` options have been removed (run `checkm qa` on the output directory for other output formats).
- New `--rm_tmp_files` option for `lineage_wf` and `taxonomy_wf` removes the `bins` and `storage` directories once the quality report has been written (other commands, such as `qa` and the plots, need these directories).

## 1.2.10 (2026-11-09)

This release includes the changes from 1.2.8 and 1.2.9, which were not released separately. **Completeness and contamination estimates can change compared with 1.2.7, but changes are expected to be rare and small.** On 1,000 randomly selected GlobDB MAGs, only 6 genomes (0.6%) changed: completeness by +0.03 for one genome and contamination by between -0.18 and +0.10 for five. No marker lineage or MIMAG quality tier changed. These changes come from two 1.2.8 bug fixes (see [#423](https://github.com/Ecogenomics/CheckM/pull/423)). The 1.2.9 and 1.2.10 changes do not affect results.

### 1.2.10

- Prodigal's error output is now captured and included in the error message when Prodigal fails (e.g., `exit code 18: Error: no input sequences to analyze.`). Results are unchanged.

### 1.2.9 (not released separately)

Performance improvements; results are unchanged.

- Marker genes are identified more efficiently: the marker file is parsed once rather than once per bin, the HMM file for each unique marker set is created once and is no longer indexed, larger bins are processed first, and HMMER runs in serial mode when bins are searched in parallel.
- Marker sets in `lineage.ms` and the reference data are parsed without `eval()`, which is ~13x faster.
- The Pfam clan file and other reference files are read once rather than once per bin, and `lineage_set` no longer searches the reference tree for each bin.
- HMMER results are parsed in parallel (`qa` now uses `-t`).
- `--pplacer_threads` now defaults to the number of threads (`-t`), to a maximum of 64. pplacer's memory usage increases only slightly with additional threads.
- README badges are cached for a shorter time.

### 1.2.8 (not released separately)

Bug fixes.

- Fix to the Pfam same-clan filter, which treated Pfams without a clan as members of the same clan. This can slightly change completeness and contamination estimates.
- Fix to reading FASTA files without a trailing newline, which dropped the last base.
- Fix to end coordinates of merged adjacent marker gene hits.
- Hits to the same marker on a run of consecutive genes (e.g., a gene split into several fragments) are now combined into a single hit, rather than being paired in an order that depended on hits to other markers. This can slightly change contamination and strain heterogeneity estimates.
- Failures of external programs (Prodigal, HMMER, pplacer) and of worker processes now cause CheckM to exit with an error instead of exiting with status 0 or silently dropping bins.
- Bins without any genes are now skipped with a warning, and FASTA records without sequence data are ignored.
- Output is now deterministic (order of HMM models, bins and markers in output files).
- Fix to infinite loop when using `--reduced_tree` for bins placed at one node of the reduced tree.
- Fix to `qa -o 4 --individual_markers` corrupting stored results, `qa -o 9` omitting markers, and parsing of custom HMM files with models lacking accessions or thresholds.
- Fix to `--ali_top_hits` aligning the worst rather than best hit.
- Fix to `gc_bias_plot` (pysam compatibility), `nx_plot` step sizes, histogram bins and percentages, and `len_hist` bar alignment.
- Fix to `ssu_finder` strand handling, contig matching, sequence coordinates and orientation of reverse-strand hits.
- Fix to handling of gzipped FASTA files in `ssu_finder`, `unique` and `modify`.
- Fix to `coverage` with BAM files lacking NM tags, BAM files with the same name, and BAM files without mapped reads.
- Fix to `checkm --help`, invalid `--tmpdir` values, bin input tables with blank lines, and `checkm data setRoot` errors.
- Paths containing spaces are now supported.
- Removed unused code (including `util/img.py` and `uniqueMarkers.py`) and fixed spelling in comments and messages.

## 1.2.7 (2026-11-08)

This release includes the changes from 1.2.6, which was not released separately.

- Fix to crash in `BinStatistics.sequenceStats` when no called genes are available.
- Fix to report hmmalign failures instead of ignoring them, which caused a `NameError` on `mask` ([#387](https://github.com/Ecogenomics/CheckM/issues/387)). *(1.2.6)*
- Fix to `OSError: AF_UNIX path too long` when using `--tmpdir` ([#402](https://github.com/Ecogenomics/CheckM/issues/402)). *(1.2.6)*
- Fix to the `coverage` command dropping bins that contain sequences with the same ID as sequences in other bins ([#386](https://github.com/Ecogenomics/CheckM/issues/386)). *(1.2.6)*

## 1.2.5 (2026-03-02)

- Fix to regular expression syntax warnings in `prettytable.py`.
- Fix to use of setuptools in `checkmData.py` (thanks to meereeum).
- Fix to pysam import in `coverage.py` (thanks to baptwr).

## 1.2.4 (2025-06-12)

- Fixed issue with Prodigal not trying the metagenomic model (`-p meta`) for failed genomes (fix provided by noelle-haddad).
- Added Dockerfile provided by Louis-MG.

## 1.2.3 (2024-06-25)

- Changes to support numpy 2.0.
- Better handling of non-zero return codes from external programs.

## 1.2.2 (2022-10-31)

- Support gzipped protein files as input (thanks to alienzj for the PR).

## 1.2.1 (2022-08-03)

- Renamed private classes to resolve a multiprocessing `AttributeError` (thanks to misialq for the PR).

## 1.2.0 (2022-04-30)

- Modified how bin IDs are identified to improve output tables when processing compressed files (e.g., the output bin ID will be `my_bin` instead of `my_bin.fna`).

## 1.1.11 (2022-04-24)

- Fixed error with `inverse_transformed` being deprecated in newer versions of Matplotlib.
- Updated minimum required versions of Python packages.

## 1.1.10 (2022-04-21)

- Fixed bug with missing `bCalledGenes` flag that was impacting a number of commands.

## 1.1.9 (2022-04-14)

- Fixed support for gzipped input FASTA files.

## 1.1.8 (2022-04-10)

- Updated `setup.py` to include a long description and fix the PyPI URL.

## 1.1.7 (2022-04-10)

- Added Jie Zhu's (alienzj) PR to allow input files to be specified via a file.
- Added Jeremy Jacobson's (jjacobson95) PR to fix a prettytable import issue.
- Added fix by Jie Li (jili6t) regarding more stringent matching of bin extensions.

## 1.1.6 (2022-04-09)

- Allow the CheckM data path to be set through the `CHECKM_DATA_PATH` environment variable.

## 1.1.5 (2022-04-09)

- Small fix to maxibor's PR to allow the path to the `DATA_CONFIG` file to be set.

## 1.1.4 (2022-04-09)

Released on PyPI only.

- Added Maxime Borry's (maxibor) PR to allow the path to the `DATA_CONFIG` file to be set.

## 1.1.3 (2020-07-09)

- Fixed error resulting from incorrect GFF files produced by Prodigal on FASTA files with Windows-style line endings.
- Replaced usage of the `normed` keyword with `density` as required for Matplotlib >= 3.2.0.

## 1.1.2 (2019-12-22)

- Removed erroneous error message when setting the root of the CheckM reference data.
- Fixed infinite loop when trying to install the CheckM reference data to a directory without write permission.

## 1.1.1 (2019-12-04)

- Fixes a few legacy Python 2 syntax issues (thanks to Sebastian Jaenicke).

## 1.1.0 (2019-12-01)

- Migration to Python 3. Massive thanks to baudrly, vinisalazar and asafpr for the initial Python 2 to 3 port.
- Fixed bug in `ssu_finder` that on rare occasions caused an SSU sequence to be reported more than once.
- Now using the i-Evalue to break ties in the rare cases when hits to a target sequence have exactly the same full E-value. This makes the code more deterministic.
- Improvements to information reported in `marker_plot`.
- Fixed CLI issue with the `outliers` command.

## 1.0.18 (2019-08-11)

- Fix to the command-line interface of `bin_qa_plot` (thanks to Finesim97 for the PR).

## 1.0.17 (2019-08-10)

- Fixed bug with the CLI of plotting commands introduced in 1.0.16.

## 1.0.16 (2019-08-08)

- Made `~/.checkm` the default directory for the CheckM reference data.
- `checkm data setRoot` is no longer interactive and the directory must be given on the command line.

## 1.0.15 (2019-08-08)

- Fixed the `checkm data setRoot` command to work from the command line (thanks to cbt for the PR).
- Added logging to files.
- Improved reporting of CheckM results to include a timestamp.
- Removed the `--force-overwrite` option; output directories are now allowed to already exist. (1.0.14, which was not released separately, had fixed `checkm test` missing this option.)

## 1.0.13 (2019-01-07)

- Failed commands now exit with error code 1.
- Added `--force-overwrite` option.

## 1.0.12 (2018-09-14)

- Fixed issue with a missing `bVerbose` parameter causing CheckM to fail on some systems; the `bVerbose` option was removed from all commands.

## 1.0.11 (2018-01-03)

- Fixed `lineage_wf` CLI bug causing a crash due to a missing `coverage_file` argument.

## 1.0.10 (2017-12-22)

- Version number update required for upload to PyPI; otherwise the same as 1.0.9.

## 1.0.9 (2017-12-22)

Released on GitHub only (see 1.0.10).

- Removed ScreamingBackpack as a dependency and modified CheckM to never remove user files.
- Added `life` as an acceptable rank for `taxonomy_wf` so the universal prokaryotic marker set can be used.
- Fixed bug with processing empty gene files.
- Fixed bug with processing identical user genomes, which caused pplacer to create polytomies.

## 1.0.8 (2017-11-09)

- Fix for the `lineage_wf` command that was causing an invalid marker set to be selected in rare cases when using `--reduced_tree`. This caused CheckM to fail with `[Error] Failed to set a selected lineage-specific marker set.`

## 1.0.7 (2016-09-02)

- Fixed incorrect reporting of the number of marker genes/sets for taxon-specific marker sets.

## 1.0.6 (2016-04-23)

- Minor bug fixes.
- Fixed `--pplacer_threads` issue with the `checkm test` command.
- Fixed issue with `checkm setRoot /my/data/directory` interactively requesting the data directory.

## 1.0.5 (2015-11-03)

- Added `--pplacer_threads` option to the `tree` and `lineage_wf` commands to allow control over the number of threads used by pplacer. This can substantially reduce the peak memory requirements of CheckM.
- Minor bug fixes.

## 1.0.4 (2015-09-11)

Minor bug fixes and enhancements.

- Removed import of pylab as this may have a GTK dependency.
- Moved to pysam >= 0.8.0 for coverage calculation.
- Added checks for expected sequence types (nucleotide vs. amino acid).
- Fixed bug with the `tetra` command when processing lowercase bases.

## 1.0.3 (2015-06-29)

- Fixed incompatibility with the GFF output of Prodigal >= 2.6.1.

## 1.0.2 (2015-06-23)

- Fixed issue with pysam that was causing the `checkm coverage` command to report no mapped reads.

## 1.0.1 (2015-06-12)

- Fixed bug in the E. coli test (`checkm test`).
- Fixed bug in the installation of CheckM from PyPI.

## 1.0.0 (2015-06-11)

- Added output format to the `qa` command that writes all marker genes to a FASTA file (courtesy of Hunter Cameron, [#31](https://github.com/Ecogenomics/CheckM/pull/31)).
- Added unique node IDs to `qa` output tables.
- Added `--exclude_markers` option to the `qa` command.
- Fix to remove temporary `*.ssi` files.

## 0.9.7 (2015-01-16)

- Added unique IDs of selected marker sets to CheckM output tables.
- Improved quality of internal taxonomic labels in the reference genome tree.

## 0.9.6 (2014-12-02)

- Fixed bug preventing `checkm test` from running.
- Moving towards PEP 8 compliance.

## 0.9.5 (2014-11-20)

- Added support for a reduced reference tree through the `--reduced_tree` option.
- Added support for using previously called genes instead of identifying these de novo.
- Improved handling of data dependencies and server connection through ScreamingBackpack.
- Improved plots for compatibility with Matplotlib 1.4.0.
- Improved command-line help.

## 0.9.4 (2014-09-09)

- Updated plots.
- Made marker set refinement on by default.
- Removed all references to local directories.
- Updated E. coli test.

## 0.9.2 (2014-09-06)

Released on PyPI only. No release notes were published.

## 0.9.1 (2014-09-05)

Released on PyPI only. No release notes were published.

## 0.9.0 (2014-09-04)

First public release. No release notes were published.
