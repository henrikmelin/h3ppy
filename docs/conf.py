# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import datetime
import os
import sys

import sphinx_rtd_theme

#sys.path.insert(0, os.path.join(os.path.split(__file__)[0], '..'))
#sys.path.insert(0, os.path.abspath('../h3ppy/'))
sys.path.insert(0, os.path.abspath('../'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'h3ppy'
copyright = f'2020–{datetime.datetime.now():%Y}, Henrik Melin'
author = 'Henrik Melin'
version = '0.4.3'
release = '0.4.3'
language = 'en'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx_rtd_theme',
    'sphinx.ext.mathjax',
]

mathjax_path = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML"


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

default_role = 'code'
master_doc = 'index'
source_suffix = [".rst", ".md"]

intersphinx_mapping = {
    'numpy': ('https://numpy.org/doc/stable', None),
}

# Autodoc
autodoc_member_order = 'bysource'
# autoclass_content = 'both'
# autodoc_typehints = 'both'
# autodoc_typehints_description_target = 'documented_params'
autodoc_inherit_docstrings = False

pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = []
# html_extra_path = ['google065a4d650d8ee82d.html']
html_logo = 'h3ppy-logo.png'
html_theme_options = {
    'logo_only': True,
}
