# Contributing to CheckM

Thank you for your interest in CheckM. Please read this page before you open an issue or a pull request.

## Project status

CheckM v1 (this repository) is **unmaintained**. Its successor is [CheckM2](https://github.com/chklovski/CheckM2).

We only make small updates to keep CheckM working with new versions of its external dependencies. These include Python, numpy, scipy, matplotlib, pysam, dendropy, setuptools, HMMER, Prodigal and gPlace. As a result:

- **We welcome** bug reports and fixes for crashes or incorrect results caused by new versions of these dependencies.
- **We will not accept** new features, new commands, changes to the marker sets or reference data, or changes to how completeness and contamination are estimated.
- Responses to issues and pull requests may be slow.

## Reporting a bug

Please report bugs through [GitHub issues](https://github.com/Ecogenomics/CheckM/issues). First check that the bug has not already been reported, and that it is present in the [latest release](https://pypi.org/project/checkm-genome/). The documentation at <https://ecogenomics.github.io/CheckM> answers many usage questions.

A useful bug report includes:

- the CheckM version (`checkm -h` prints it) and how CheckM was installed (pip, Bioconda, Docker)
- your Python version and operating system
- the versions of HMMER, Prodigal and gPlace (or pplacer, for CheckM < 1.3.0)
- the exact command you ran
- the full error message, along with `checkm.log` from the output directory
- if possible, a small set of genomes that reproduces the problem

## Setting up a development environment

CheckM requires Python >= 3.11. You also need the following on your `PATH`:

- [HMMER](http://hmmer.org/) >= 3.1
- [Prodigal](https://github.com/hyattpd/Prodigal) 2.6.3
- [gPlace](https://github.com/cmc-aau/gPlace) >= 0.1.0 (`pip install gtdb-gplace`; requires Python >= 3.12 and can be installed in a separate environment)

Install CheckM from your clone in editable mode:

```bash
git clone https://github.com/Ecogenomics/CheckM.git
cd CheckM
pip install -e .
```

You also need the CheckM reference data (`checkm_data_2015_01_16.tar.gz`). See [Reference data](https://ecogenomics.github.io/CheckM/installing/reference_data.html).

**Set `CHECKM_DATA_PATH` before running CheckM or its tests.** If it is unset, importing `checkm` reads the data location from `checkm/DATA_CONFIG`. If no data is found there, CheckM prompts for a location and writes it into `checkm/DATA_CONFIG`. Never commit that change: `dataRoot` must remain `""` in the repository.

```bash
export CHECKM_DATA_PATH=/path/to/checkm_data
```

## Running the tests

```bash
python3 -m unittest discover -s checkm/test -t .                 # all unit tests
python3 -m unittest checkm.test.test_markerSets                  # one module
python3 -m unittest checkm.test.test_markerSets.VerifyMarkerSet  # one test class
checkm test <output_dir>                                         # end-to-end run on an E. coli genome
```

Tests that need an external program are skipped when it isn't on your `PATH`. Before you submit a pull request, make sure that the unit tests and `checkm test` both pass.

## Making changes

- Keep each pull request small and focused on a single problem.
- If a fix is not small, open an issue to discuss it first.
- Follow the style of the code around your change. The project has no formatter or linter configuration.
- Add a unit test under `checkm/test/` that fails without your fix.
- Results should not change. If a fix changes completeness, contamination or other output, explain why in the pull request.
- If you change a command's arguments or help text, update both `bin/checkm` (the argument parser and the help text in `printHelp()`) and `checkm/main.py`. Then regenerate the documentation's help files as described in [docs/README.md](docs/README.md).
- Commit messages use a short type prefix, such as `fix: ...`, `perf: ...`, `docs: ...` or `chore: ...`.

## Pull requests

Open pull requests against the `master` branch. Please include:

- the problem being fixed, with a link to the issue if there is one
- the dependency versions involved
- how you tested the change

Do not change the version number. The maintainers set it when they make a release.

## Releases (maintainers)

Each fix is released as a point release:

1. Set the new version number in `checkm/VERSION`, which holds only the version number. Add a section for the release at the top of [CHANGELOG.md](CHANGELOG.md), which is also shown on the documentation site.
2. Merge to `master`. This rebuilds the [documentation site](https://ecogenomics.github.io/CheckM).
3. Tag the release as `vX.Y.Z` and create a GitHub release.
4. Build and upload the package to [PyPI](https://pypi.org/project/checkm-genome/) as `checkm-genome`. The Bioconda recipe is updated separately.

## License

CheckM is licensed under the GNU General Public License v3 (see [LICENSE](LICENSE)). By contributing, you agree that your contributions will be licensed under the same terms.

## Acknowledgements

CheckM was developed by Donovan Parks, Michael Imelfort and Connor Skennerton, with Philip Hugenholtz and Gene Tyson, at the Australian Centre for Ecogenomics, University of Queensland. If you use CheckM, please [cite it](https://ecogenomics.github.io/CheckM/references.html).

### Python 3

Many thanks to [Vini Salazar](https://github.com/vinisalazar), the primary developer of the CheckM Python 3 branch. Thanks also to [baudrly](https://github.com/baudrly) and [Asaf Peer](https://github.com/asafpr), who contributed to the initial port from Python 2 to Python 3, and to Sebastian Jaenicke, who fixed remaining Python 2 syntax.

### Code contributions

We also thank the following people for contributing code:

- Ben Woodcroft: the `bin_union` command
- Hunter Cameron: output format 10 for `qa`, which writes a FASTA file of the marker genes in each bin
- alienzj: gzipped protein files and tables of bins as input
- Noëlle Haddad: a fix so that Prodigal's metagenomic mode is tried for genomes that fail
- meereeum: compatibility with old and new versions of setuptools
- baptwr: a fix to the pysam import used by the coverage commands
- Louis-MG: the Dockerfile
- C. Titus Brown and Maxime Borry: fixes to `checkm data setRoot`
- Jie Li: stricter matching of bin file extensions
- Finesim97: a fix to the `bin_qa_plot` command
- Matthias Bernt: names for all subcommand parsers
- Michal Ziemski: renaming of private class methods
- Jeremy Jacobson: a fix to the vendored PrettyTable

Other contributions are acknowledged in the [change log](CHANGELOG.md).

### Community

We are grateful to the community members who package and maintain CheckM on [Bioconda](https://anaconda.org/bioconda/checkm-genome). We also thank everyone who has reported bugs and suggested fixes through GitHub issues.

CheckM relies on [HMMER](http://hmmer.org/), [Prodigal](https://github.com/hyattpd/Prodigal), [pplacer](https://matsen.fhcrc.org/pplacer/) and [gPlace](https://github.com/cmc-aau/gPlace), and we thank their developers.
