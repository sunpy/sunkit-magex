# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

import os

# -- Project information -----------------------------------------------------

# The full version, including alpha/beta/rc tags
from sunkit_magex import __version__
release = __version__

project = "sunkit-magex"
copyright = "2024, The SunPy Community"
author = "The SunPy Community"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named "sphinx.ext.*") or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
    "sphinx_automodapi.automodapi",
    "sphinx_automodapi.smart_resolver",
    'sphinx_gallery.gen_gallery',
]

# Add any paths that contain templates here, relative to this directory.
# templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# -- Options for intersphinx extension ---------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/", None),
    "astropy": ("https://docs.astropy.org/en/stable", None),
    "sunpy": ("https://docs.sunpy.org/en/stable", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "streamtracer": ("https://streamtracer.readthedocs.io/en/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "reproject": ("https://reproject.readthedocs.io/en/stable/", None),
}

# -- Options for HTML output -------------------------------------------------

html_theme = "sunpy"

# -- Sphinx Gallery ----------------------------------------------------------

from sphinx_gallery.sorting import ExplicitOrder  # noqa

sphinx_gallery_conf = {
    "ignore_pattern": ".*helpers.py",
    "examples_dirs": "../examples",
    "gallery_dirs": os.path.join("generated", "gallery"),
    "subsection_order": ExplicitOrder(["../examples/using_pfsspy",
                                       "../examples/finding_data",
                                       "../examples/utils",
                                       "../examples/pfsspy_info",
                                       "../examples/testing"]),
    "reference_url": {"sphinx_gallery": None}
}

# -- Other options ----------------------------------------------------------
default_role = 'py:obj'

os.environ["JSOC_EMAIL"] = 'jsoc@sunpy.org'

nitpicky = True
