# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

import os
import datetime

from packaging.version import Version

# -- Project information -----------------------------------------------------

# The full version, including alpha/beta/rc tags
from sunkit_magex import __version__

_version = Version(__version__)
version = release = str(_version)
# Avoid "post" appearing in version string in rendered docs
if _version.is_postrelease:
    version = release = _version.base_version
# Avoid long githashes in rendered Sphinx docs
elif _version.is_devrelease:
    version = release = f"{_version.base_version}.dev{_version.dev}"
is_development = _version.is_devrelease
is_release = not(_version.is_prerelease or _version.is_devrelease)

project = "sunkit-magex"
author = "The SunPy Community"
copyright = f"{datetime.datetime.now().year}, {author}"  # noqa: A001

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

# Treat everything in single ` as a Python reference.
default_role = "py:obj"

# -- Options for intersphinx extension ---------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/", None),
    "astropy": ("https://docs.astropy.org/en/stable", None),
    "sunpy": ("https://docs.sunpy.org/en/stable", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "streamtracer": ("https://streamtracer.readthedocs.io/en/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "reproject": ("https://reproject.readthedocs.io/en/stable/", None),
    "pfsspy": ("https://pfsspy.readthedocs.io/en/latest/", None),
}

# -- Options for HTML output -------------------------------------------------

html_theme = "sunpy"

# Render inheritance diagrams in SVG
graphviz_output_format = "svg"

graphviz_dot_args = [
    "-Nfontsize=10",
    "-Nfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Efontsize=10",
    "-Efontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Gfontsize=10",
    "-Gfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]  # NOQA: ERA001

# By default, when rendering docstrings for classes, sphinx.ext.autodoc will
# make docs with the class-level docstring and the class-method docstrings,
# but not the __init__ docstring, which often contains the parameters to
# class constructors across the scientific Python ecosystem. The option below
# will append the __init__ docstring to the class-level docstring when rendering
# the docs. For more options, see:
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#confval-autoclass_content
autoclass_content = "both"

# -- Other options ----------------------------------------------------------

# JSOC email OS ENV
# See https://github.com/sunpy/sunpy/wiki/Home:-JSOC
os.environ["JSOC_EMAIL"] = 'jsoc@sunpy.org'

nitpicky = True
numfig = True

# -- Sphinx Gallery ----------------------------------------------------------

from sunpy_sphinx_theme import PNG_ICON  # noqa
from sphinx_gallery.sorting import ExplicitOrder  # noqa

sphinx_gallery_conf = {
    "abort_on_example_error": False,
    "backreferences_dir": os.path.join("generated", "modules"),
    "default_thumb_file": PNG_ICON,
    "examples_dirs": os.path.join("..", "examples"),
    "filename_pattern": '^((?!skip_).)*$',
    "gallery_dirs": os.path.join("generated", "gallery"),
    "ignore_pattern": "helpers.py",
    "matplotlib_animations": True,
    "only_warn_on_example_error": True,
    "plot_gallery": True,
    "remove_config_comments": True,
    "subsection_order": ExplicitOrder([
        "../examples/using_sunkit_magex_pfss",
        "../examples/finding_data",
        "../examples/utils",
        "../examples/internals",
        "../examples/testing"
    ]),
}
