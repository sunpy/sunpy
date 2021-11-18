"""
Configuration file for the Sphinx documentation builder.
"""
# -- stdlib imports ------------------------------------------------------------
import os
import sys
import datetime
import warnings
from pkg_resources import get_distribution
from packaging.version import Version

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
    os.environ['HIDE_PARFIVE_PROGESS'] = 'True'

# -- Non stdlib imports --------------------------------------------------------
import ruamel.yaml as yaml  # NOQA
from sphinx_gallery.sorting import ExplicitOrder  # NOQA
from sphinx_gallery.sorting import ExampleTitleSortKey  # NOQA

import sunpy  # NOQA
from sunpy import __version__  # NOQA
from sunpy.util.exceptions import SunpyDeprecationWarning, SunpyPendingDeprecationWarning  # NOQA
from matplotlib import MatplotlibDeprecationWarning  # NOQA
from astropy.utils.exceptions import AstropyDeprecationWarning  # NOQA
# -- Project information -------------------------------------------------------
project = 'SunPy'
author = 'The SunPy Community'
copyright = '{}, {}'.format(datetime.datetime.now().year, author)


# Register remote data option with doctest
import doctest  # NOQA
REMOTE_DATA = doctest.register_optionflag('REMOTE_DATA')

# The full version, including alpha/beta/rc tags
release = __version__
sunpy_version = Version(__version__)
is_release = not(sunpy_version.is_prerelease or sunpy_version.is_devrelease)

# We want to ignore all warnings in a release version.
if is_release:
    warnings.simplefilter("ignore")
warnings.filterwarnings("error", category=SunpyDeprecationWarning)
warnings.filterwarnings("error", category=SunpyPendingDeprecationWarning)
warnings.filterwarnings("error", category=MatplotlibDeprecationWarning)
warnings.filterwarnings("error", category=AstropyDeprecationWarning)

# -- SunPy Sample Data and Config ----------------------------------------------
# We set the logger to debug so that we can see any sample data download errors
# in the CI, especially RTD.
ori_level = sunpy.log.level
sunpy.log.setLevel("DEBUG")
import sunpy.data.sample  # NOQA
sunpy.log.setLevel(ori_level)

# For the linkcheck
linkcheck_ignore = [r"https://doi.org/\d+",
                    r"https://element.io/\d+",
                    r"https://github.com/\d+",
                    r"https://docs.sunpy.org/\d+"]
linkcheck_anchors = False

# -- General configuration -----------------------------------------------------
# sphinxext-opengraph
ogp_image = "https://raw.githubusercontent.com/sunpy/sunpy-logo/master/generated/sunpy_logo_word.png"
ogp_use_first_image = True
ogp_description_length = 160
ogp_custom_meta_tags = [
    '<meta property="og:ignore_canonical" content="true" />',
]

# Suppress warnings about overriding directives as we overload some of the
# doctest extensions.
suppress_warnings = ['app.add_directive', ]

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'matplotlib.sphinxext.plot_directive',
    'sphinx_automodapi.automodapi',
    'sphinx_automodapi.smart_resolver',
    'sphinx_changelog',
    'sphinx_gallery.gen_gallery',
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.inheritance_diagram',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sunpy.util.sphinx.doctest',
    'sunpy.util.sphinx.generate',
    "sphinxext.opengraph",
    'sphinx_panels',
]

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

if is_release:
    exclude_patterns.append('dev_guide/contents/*')

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

# Enable nitpicky mode, which forces links to be non-broken
nitpicky = True
# This is not used. See docs/nitpick-exceptions file for the actual listing.
nitpick_ignore = []
for line in open('nitpick-exceptions'):
    if line.strip() == "" or line.startswith("#"):
        continue
    dtype, target = line.split(None, 1)
    target = target.strip()
    nitpick_ignore.append((dtype, target))


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
    "aiapy": ("https://aiapy.readthedocs.io/en/stable/", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "astroquery": ("https://astroquery.readthedocs.io/en/latest/", None),
    "drms": ("https://docs.sunpy.org/projects/drms/en/stable/", None),
    "mpl_animators": ("https://docs.sunpy.org/projects/mpl-animators/en/stable/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "parfive": ("https://parfive.readthedocs.io/en/stable/", None),
    "reproject": ("https://reproject.readthedocs.io/en/stable/", None),
    "skimage": ("https://scikit-image.org/docs/stable/", None),
    "sqlalchemy": ("https://docs.sqlalchemy.org/en/latest/", None),
    "sunkit_instruments": ("https://docs.sunpy.org/projects/sunkit-instruments/en/stable/", None),
    "zeep": ("https://docs.python-zeep.org/en/stable/", None),
}

# -- Options for HTML output ---------------------------------------------------
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

from sunpy_sphinx_theme.conf import *  # NOQA

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
# JSOC email os env
os.environ["JSOC_EMAIL"] = "jsoc@cadair.com"
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
    'plot_gallery': 'True',
    'remove_config_comments': True,
    'doc_module': ('sunpy'),
    'only_warn_on_example_error': True,
}

# -- Stability Page ------------------------------------------------------------
with open('./code_ref/sunpy_stability.yaml', 'r') as estability:
    sunpy_modules = yaml.load(estability.read(), Loader=yaml.Loader)

html_context = {
    'sunpy_modules': sunpy_modules,
    'is_development': not is_release,
}


def rstjinja(app, docname, source):
    """
    Render our pages as a jinja template for fancy templating goodness.
    """
    # Make sure we're outputting HTML
    if app.builder.format != 'html':
        return
    files_to_render = ["code_ref/stability", "dev_guide/index"]
    if docname in files_to_render:
        print(f"Jinja rendering {docname}")
        rendered = app.builder.templates.render_string(
            source[0], app.config.html_context
        )
        source[0] = rendered


# -- Sphinx setup --------------------------------------------------------------
def setup(app):
    # Generate the stability page
    app.connect("source-read", rstjinja)
