# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
# SunPy documentation build configuration file.
#
# This file is execfile()d with the current directory set to its containing
# dir.
#
# Note that not all possible configuration values are present in this file.
#
# All configuration values have a default. Some values are defined in
# the global Astropy configuration which is loaded here before anything else.
# See astropy.sphinx.conf for which values are set there.

import os
import sys
import datetime

# -- Load astropy_helpers -----------------------------------------------------

try:
    # Has astropy_helpers been installed via pip or similar?
    import astropy_helpers
except ImportError:
    # Building from the doc/source directory?
    if os.path.basename(os.getcwd()) == 'source':
        a_h_path = os.path.abspath(os.path.join('..', '..', 'astropy_helpers'))
        if os.path.isdir(a_h_path):
            sys.path.insert(1, a_h_path)

# -- Read the Docs Setup  -----------------------------------------------------

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if on_rtd:
    os.environ['SUNPY_CONFIGDIR'] = '/home/docs/'
    os.environ['HOME'] = '/home/docs/'

# -- Download Sample Data -----------------------------------------------------

import sunpy.data
sunpy.data.download_sample_data(overwrite=False)

# -- General configuration ----------------------------------------------------


# Load all of the global Astropy configuration
from astropy_helpers.sphinx.conf import *

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = '1.1'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns.append('_templates')

# This is added to the end of RST files - a good place to put substitutions to
# be used globally.
rst_epilog += """
"""

# put auto generated api docs in the generated folder.
automodapi_toctreedirnm = 'generated/api'

# -- Project information ------------------------------------------------------

# This does not *have* to match the package name, but typically does
project = u'SunPy'
author = u'The SunPy Community'
copyright = u'{}, {}'.format(datetime.datetime.now().year, author)

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.

import sunpy

# The short X.Y version.
version = sunpy.__version__.split('-', 1)[0]
# The full version, including alpha/beta/rc tags.
release = sunpy.__version__

intersphinx_mapping.pop('h5py', None)
intersphinx_mapping['astropy'] = ('http://docs.astropy.org/en/stable/', None)
intersphinx_mapping['sqlalchemy'] = ('http://docs2.sqlalchemy.org/en/latest/', None)
intersphinx_mapping['pandas'] = ('http://pandas.pydata.org/pandas-docs/stable/', None)
intersphinx_mapping['skimage'] = ('http://scikit-image.org/docs/stable/', None)
intersphinx_mapping['wcsaxes'] = ('http://wcsaxes.readthedocs.io/en/stable/', None)

# -- Options for HTML output ---------------------------------------------------

try:
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
except ImportError:
    html_theme = 'default'

# Custom sidebar templates, maps document names to template names.
# html_sidebars = {}

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "../logo/favicon.ico"
# html_logo = "../logo/sunpy_logo_compact_192x239.png"

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
# html_last_updated_fmt = ''

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = '{0} v{1}'.format(project, release)

# Output file base name for HTML help builder.
htmlhelp_basename = project + 'doc'


# -- Options for LaTeX output -------------------------------------------------

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [('index', project + '.tex', project + u' Documentation',
                    author, 'manual')]


# -- Options for manual page output -------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [('index', project.lower(), project + u' Documentation',
              [author], 1)]

# -- Swap to Napoleon ---------------------------------------------------------

extensions.remove('astropy_helpers.sphinx.ext.numpydoc')
extensions.append('sphinx.ext.napoleon')

# Disable having a separate return type row
napoleon_use_rtype = False
# Disable google style docstrings
napoleon_google_docstring = False


# -- Options for the edit_on_github extension ---------------------------------
extensions += ['astropy_helpers.sphinx.ext.edit_on_github',
               'sphinx.ext.doctest']


# Don't import the module as "version" or it will override the
# "version" configuration parameter
# TODO: make this smart like astropy
edit_on_github_project = "sunpy/sunpy"
edit_on_github_branch = "master"

edit_on_github_source_root = ""
edit_on_github_doc_root = "docs"


# -- Sphinx Gallery ----------------------------------------------------------
try:
    import sphinx_gallery
    extensions += ['sphinx_gallery.gen_gallery']

    sphinx_gallery_conf = {
        # path to store the module using example template
        'mod_example_dir': 'generated{}modules'.format(os.sep),
        # execute all examples except those that start with "skip_"
        'filename_pattern': '^((?!skip_).)*$',
        # path to the examples scripts
        'examples_dirs': '..{}..{}examples'.format(os.sep, os.sep),
        'gallery_dirs': 'generated{}gallery'.format(os.sep),
        'default_thumb_file': '..{}logo{}sunpy_icon_128x128.png'.format(os.sep, os.sep),
        'reference_url': {
            'sunpy': None,
            'astropy': 'http://docs.astropy.org/en/stable/',
            'matplotlib': 'http://matplotlib.org/',
            'numpy': 'http://docs.scipy.org/doc/numpy/'}
        }
except ImportError:
    def setup(app):
        app.warn('The sphinx_gallery extension is not installed, so the '
                 'gallery will not be built.  You will probably see '
                 'additional warnings about undefined references due '
                 'to this.')
