# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import sphinx_drove_theme 

sys.path.insert(0, os.path.abspath('../../src/PyGran'))


# -- Project information -----------------------------------------------------

project = 'PyGran'
copyright = '2019, Andrew Abi-Mansour'
author = 'Andrew Abi-Mansour'

# The full version, including alpha/beta/rc tags
release = '1.2.0'

# -- General configuration ------------------------------------------------
autoclass_content = "both"  # include both class docstring and __init__
autodoc_member_order = 'bysource'
autodoc_default_flags = [
	# Make sure that any autodoc declarations show the right members
	"members",
	"inherited-members",
	"private-members",
	"show-inheritance"]

autosummary_generate = True  # Make _autosummary files and include them
napoleon_numpy_docstring = False  # Force consistency, leave only Google
napoleon_use_rtype = False  # More legible
todo_include_todos = True

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 
	'sphinx.ext.autosummary',
	'sphinx.ext.viewcode',
	# The Napoleon extension allows for nicer argument formatting.
	'sphinx.ext.napoleon',
	'sphinxcontrib.bibtex',
	'matplotlib.sphinxext.plot_directive',
	'sphinx.ext.imgmath',
	'sphinx.ext.todo']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_drove_theme' 
html_theme_path = [sphinx_drove_theme.get_html_theme_path()]

numfig = True
math_number_all = True
