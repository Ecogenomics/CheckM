# CheckM documentation

Source for the CheckM documentation site, https://ecogenomics.github.io/CheckM/.
It is built with [Sphinx](https://www.sphinx-doc.org/) and the Read the Docs theme.
The `.github/workflows/docs.yml` workflow builds the site and publishes it to the
`gh-pages` branch when files under `docs/` change on `master`.

## Building locally

```bash
cd docs
python -m pip install -r requirements.txt
make html SPHINXOPTS="-W --keep-going"   # warnings are errors, as in the workflow
```

Open `build/html/index.html` in a browser. The version shown is read from
`../checkm/VERSION`, and the change log page includes `../CHANGELOG.md`
(rendered with MyST).

## Regenerating the command help

The `Arguments` section of each command page includes the output of
`checkm <command> -h` from `src/commands/help/<command>.txt`
(`src/commands/help/checkm.txt` holds `checkm -h`). Regenerate these files
after changing a command's arguments or help text:

```bash
export CHECKM_DATA_PATH=/path/to/checkm_data   # required; see update_help.sh
./update_help.sh /path/to/checkm               # defaults to checkm on PATH
```

Run it with a CheckM installed from this source tree (e.g. `pip install -e .`)
so the help matches the code. The help text was last generated with Python 3.11;
the formatting of argparse help (e.g. `options:` vs `optional arguments:`)
depends on the Python version.

When adding a command, also add a page under `src/commands/` and list it in
`src/commands/index.rst`.
