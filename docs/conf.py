# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'LineageGRN'
copyright = '2025, Xiaochen Yu'
author = 'Xiaochen Yu'
release = '0.0.1'


import os
import sys
sys.path.append(os.path.relpath('./source'))
sys.path.append(os.path.relpath('.'))

extensions = [
    'myst_parser',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
]



source_suffix={
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
autodoc_member_order = 'bysource'

import os
if not os.path.exists('_static'):
    os.makedirs('_static')

html_theme = "furo"  
html_static_path = ["_static"]
html_css_files = [
    "css/override.css", 
]
html_css_files = [
    'custom.css',
]

html_logo = "_static/logo.png"

html_theme_options = {
    "sidebar_hide_name": True,
    "light_css_variables": {
        "color-brand-primary": "#357473",
        "color-brand-content": "#357473",
    },
}



html_theme = 'furo'
html_static_path = ['_static']




