import os
import sys

# The pyEDW package lives one level above docs/; its parent goes on the path so
# ``import pyEDW`` resolves for autodoc.
sys.path.insert(0, os.path.abspath("../.."))

project = "pyEDW"
author = "Damian R. Sowinski"
copyright = "2026, Damian R. Sowinski"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
]

napoleon_numpy_docstring = True
napoleon_google_docstring = False

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

autodoc_member_order = "bysource"
autodoc_mock_imports = ["numba"]
