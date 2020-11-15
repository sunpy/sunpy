"""
Configuration file for the Sphinx documentation builder.

isort:skip_file
"""
# flake8: NOQA: E402

# -- stdlib imports ------------------------------------------------------------
from sunpy_sphinx_theme.conf import *
import sunpy.data.sample  # isort:skip
from sunpy import __version__
import sunpy
from sphinx_gallery.sorting import ExampleTitleSortKey
from sphinx_gallery.sorting import ExplicitOrder
import ruamel.yaml as yaml
import os
import sys
import datetime
from pkg_resources import get_distribution, DistributionNotFound

# -- Check for dependencies ----------------------------------------------------

doc_requires = get_distribution("sunpy").requires(extras=("docs",))
missing_requirements = []
for requirement in doc_requires:
    try:
        get_distribution(requirement)
    except Exception as e:
        missing_requirements.append(requirement.name)
if missing_requirements:
    print(
        f"The {' '.join(missing_requirements)} package(s) could not be found and "
        "is needed to build the documentation, please install the 'docs' requirements."
    )
    sys.exit(1)

# -- Read the Docs Specific Configuration --------------------------------------

# This needs to be done before sunpy is imported
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd:
    os.environ['SUNPY_CONFIGDIR'] = '/home/docs/'
    os.environ['HOME'] = '/home/docs/'
    os.environ['LANG'] = 'C'
    os.environ['LC_ALL'] = 'C'

# -- Non stdlib imports --------------------------------------------------------


# -- Project information -------------------------------------------------------

project = 'SunPy'
author = 'The SunPy Community'
copyright = '{}, {}'.format(datetime.datetime.now().year, author)

# The full version, including alpha/beta/rc tags
release = __version__
is_development = '.dev' in __version__

# -- SunPy Sample Data and Config ----------------------------------------------

# We set the logger to debug so that we can see any sample data download errors
# in the CI, especially RTD.
ori_level = sunpy.log.level
sunpy.log.setLevel("DEBUG")
sunpy.log.setLevel(ori_level)

# For the linkcheck
linkcheck_ignore = [r"https://doi.org/\d+",
                    r"https://riot.im/\d+",
                    r"https://github.com/\d+",
                    r"https://docs.sunpy.org/\d+"]
linkcheck_anchors = False

# This is added to the end of RST files - a good place to put substitutions to
# be used globally.
rst_epilog = """
.. SunPy
.. _SunPy: https://sunpy.org
.. _`SunPy mailing list`: https://groups.google.com/group/sunpy
.. _`SunPy dev mailing list`: https://groups.google.com/group/sunpy-dev
"""

# -- General configuration -----------------------------------------------------

# Suppress warnings about overriding directives as we overload some of the
# doctest extensions.
suppress_warnings = ['app.add_directive', ]

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.inheritance_diagram',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.doctest',
    'sphinx.ext.mathjax',
    'sphinx_automodapi.automodapi',
    'sphinx_automodapi.smart_resolver',
    'sphinx_gallery.gen_gallery',
    'sunpy.util.sphinx.doctest',
    'matplotlib.sphinxext.plot_directive',
    'sunpy.util.sphinx.minigallery',
    'sunpy.util.sphinx.generate',
    'sunpy.util.sphinx.changelog',
]

if on_rtd:
    rtd_version = os.environ.get('READTHEDOCS_VERSION', None)
    stable = (rtd_version == 'stable') or (rtd_version[0] == 'v')
    if stable:
        # Stable versions have a manually rendered changelog, so remove the
        # changelog extension
        extensions.remove('sunpy.util.sphinx.changelog')

# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.

# Add any extra paths that contain custom files (such as robots.txt or
# .htaccess) here, relative to this directory. These files are copied
# directly to the root of the documentation.
html_extra_path = ['robots.txt']

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The reST default role (used for this markup: `text`) to use for all
# documents. Set to the "smart" one.
default_role = 'obj'

# Disable having a separate return type row
napoleon_use_rtype = False

# Disable google style docstrings
napoleon_google_docstring = False

# -- Options for intersphinx extension -----------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "python": (
        "https://docs.python.org/3/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/python3.inv"),
    ),
    "numpy": (
        "https://numpy.org/doc/stable/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/numpy.inv"),
    ),
    "scipy": (
        "https://docs.scipy.org/doc/scipy/reference/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/scipy.inv"),
    ),
    "matplotlib": (
        "https://matplotlib.org/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/matplotlib.inv"),
    ),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "sqlalchemy": ("https://docs.sqlalchemy.org/en/latest/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "skimage": ("https://scikit-image.org/docs/stable/", None),
    "drms": ("https://docs.sunpy.org/projects/drms/en/stable/", None),
    "parfive": ("https://parfive.readthedocs.io/en/latest/", None),
    "reproject": ("https://reproject.readthedocs.io/en/stable/", None),
    "aiapy": ("https://aiapy.readthedocs.io/en/stable/", None),
}

# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# Render inheritance diagrams in SVG
graphviz_output_format = "svg"

graphviz_dot_args = [
    '-Nfontsize=10',
    '-Nfontname=Helvetica Neue, Helvetica, Arial, sans-serif',
    '-Efontsize=10',
    '-Efontname=Helvetica Neue, Helvetica, Arial, sans-serif',
    '-Gfontsize=10',
    '-Gfontname=Helvetica Neue, Helvetica, Arial, sans-serif'
]

# -- Sphinx Gallery ------------------------------------------------------------

sphinx_gallery_conf = {
    'backreferences_dir': os.path.join('generated', 'modules'),
    'filename_pattern': '^((?!skip_).)*$',
    'examples_dirs': os.path.join('..', 'examples'),
    'subsection_order': ExplicitOrder([
        '../examples/acquiring_data',
        '../examples/map',
        '../examples/map_transformations',
        '../examples/time_series',
        '../examples/units_and_coordinates',
        '../examples/plotting',
        '../examples/differential_rotation',
        '../examples/saving_and_loading_data',
        '../examples/computer_vision_techniques',
        '../examples/developer_tools'
    ]),
    'within_subsection_order': ExampleTitleSortKey,
    'gallery_dirs': os.path.join('generated', 'gallery'),
    # Comes from the theme.
    "default_thumb_file": os.path.join(html_static_path[0], "img", "sunpy_icon_128x128.png"),
    'abort_on_example_error': False,
    'plot_gallery': True,
    'remove_config_comments': True,
    'doc_module': ('sunpy')
}

# -- Stability Page ------------------------------------------------------------

with open('./dev_guide/sunpy_stability.yaml', 'r') as estability:
    sunpy_modules = yaml.load(estability.read(), Loader=yaml.Loader)

html_context = {
    'sunpy_modules': sunpy_modules
}


def rstjinja(app, docname, source):
    """
    Render our pages as a jinja template for fancy templating goodness.
    """
    # Make sure we're outputting HTML
    if app.builder.format != 'html':
        return
    src = source[0]
    if "Current status" in src[:20]:
        rendered = app.builder.templates.render_string(
            src, app.config.html_context
        )
        source[0] = rendered

# -- Sphinx setup --------------------------------------------------------------


def setup(app):
    # Generate the stability page
    app.connect("source-read", rstjinja)

    # The theme conf provides a fix for circle ci redirections
    fix_circleci(app)
