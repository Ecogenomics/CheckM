# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
from datetime import datetime

# -- Project information -----------------------------------------------------

# The version is read from checkm/VERSION (first line) rather than by importing
# the checkm package, as importing checkm requires the CheckM reference data.
_version_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '..', '..', 'checkm', 'VERSION')
with open(_version_file) as f:
    _version = f.readline().strip()

project = 'CheckM'
author = 'Donovan Parks'
copyright = f'2014-{datetime.now().year}, Donovan Parks, Connor Skennerton, Michael Imelfort'

release = _version
version = _version

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- General configuration ---------------------------------------------------

# myst_parser renders CHANGELOG.md, which changelog.rst includes
extensions = ['sphinx_sitemap', 'myst_parser']

# changelog.rst includes CHANGELOG.md without its "# Change log" title, so the included
# headings start at H2 under the page title given in changelog.rst
suppress_warnings = ['myst.header']

templates_path = []
exclude_patterns = []

# copied to the root of the built site
html_extra_path = ['robots.txt']

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'logo_only': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

html_static_path = ['_static']

html_css_files = [
    'theme_overrides.css',
]

# "Edit on GitHub" link
html_context = {
    'display_github': True,
    'github_user': 'Ecogenomics',
    'github_repo': 'CheckM',
    'github_version': 'master/docs/src/'
}

html_logo = '_static/checkm_logo.svg'
html_favicon = '_static/checkm_icon.png'

# Sitemap settings
html_baseurl = 'https://ecogenomics.github.io/CheckM/'
sitemap_url_scheme = '{link}'
