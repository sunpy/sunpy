# flake8: noqa
import os
import sys
from pathlib import Path
import datetime
from docutils import nodes, statemachine
from docutils.io import FileInput
from docutils.parsers.rst import Directive, directives
from io import StringIO

# -- Import base config from sphinx-astropy ------------------------------------
try:
    from sphinx_astropy.conf.v1 import *
except ImportError:
    print('ERROR: the documentation requires the "sphinx-astropy" package to be installed')
    sys.exit(1)

try:
    import sphinx_gallery
    if on_rtd and os.environ.get('READTHEDOCS_PROJECT').lower() != 'sunpy':
        # Gallery takes too long on RTD to build unless you have extra build time.
        has_sphinx_gallery = False
    else:
        has_sphinx_gallery = True
except ImportError:
    has_sphinx_gallery = False

if on_rtd:
    os.environ['SUNPY_CONFIGDIR'] = '/home/docs/'
    os.environ['HOME'] = '/home/docs/'
    os.environ['LANG'] = 'C'
    os.environ['LC_ALL'] = 'C'

try:
    import zeep
except ImportError as e:
    raise Exception(e, 'ERROR: zeep could not be imported. Building the documentation requires '
                    'the "zeep" package to be installed')

try:
    import skimage
except ImportError as e:
    raise Exception(e, 'ERROR: skimage could not be imported. Building the documentation requires '
                    'the "scikit-image" package to be installed')

try:
    import drms
except ImportError as e:
    raise Exception(e, 'ERROR: drms could not be imported. Building the documentation requires '
                    'the "drms" package to be installed')

try:
    import glymur
except ImportError as e:
    raise Exception(e, 'ERROR: glymur could not be imported. Building the documentation requires '
                    'the "glymur" package to be installed')

try:
    import sqlalchemy
except ImportError as e:
    raise Exception(e, 'ERROR: sqlalchemy could not be imported. Building the documentation requires '
                    'the "sqlalchemy" package to be installed')

try:
    import astroquery
except ImportError as e:
    raise Exception(e, 'ERROR: astroquery could not be imported. Building the documentation requires '
                    'the "astroquery" package to be installed')

try:
    import jplephem
except ImportError as e:
    raise Exception(e, 'ERROR: jplephem could not be imported. Building the documentation requires '
                    'the "jplephem" package to be installed')

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
# The short X.Y version.
from sunpy import __version__  # noqa isort:skip
version = '.'.join(__version__.split('.')[:3])
# The full version, including alpha/beta/rc tags.
release = __version__
# Is this version a development release
is_development = '.dev' in release

# -- Shut up numpy warnings from WCSAxes --------------------------------------
import numpy as np  # noqa isort:skip
np.seterr(invalid='ignore')

# -- Download Sample Data -----------------------------------------------------
import sunpy.data.sample  # noqa isort:skip

# -- General configuration ----------------------------------------------------
# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = '2.0'

# To perform a Sphinx version check that needs to be more specific than
# major.minor, call `check_sphinx_version("x.y.z")` here.
check_sphinx_version(needs_sphinx)

# Add any custom intersphinx for SunPy
intersphinx_mapping.pop('h5py', None)
intersphinx_mapping['sqlalchemy'] = ('https://docs.sqlalchemy.org/en/latest/', None)
intersphinx_mapping['pandas'] = ('https://pandas.pydata.org/pandas-docs/stable/', None)
intersphinx_mapping['skimage'] = ('https://scikit-image.org/docs/stable/', None)
intersphinx_mapping['drms'] = ('https://docs.sunpy.org/projects/drms/en/stable/', None)
intersphinx_mapping['parfive'] = ('https://parfive.readthedocs.io/en/latest/', None)
intersphinx_mapping['reproject'] = ('https://reproject.readthedocs.io/en/stable/', None)

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns.append('_templates')

# Add any paths that contain templates here, relative to this directory.
if 'templates_path' not in locals():  # in case parent conf.py defines it
    templates_path = []
templates_path.append('_templates')

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

# -- Project information ------------------------------------------------------
project = 'SunPy'
author = 'The SunPy Community'
copyright = '{}, {}'.format(datetime.datetime.now().year, author)

fix_circleci = lambda x: None
try:
    from sunpy_sphinx_theme.conf import *
except ImportError:
    html_theme = 'default'

try:
    import ruamel.yaml as yaml
    has_yaml = True
    # Load data about stability
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
except ImportError:
    has_yaml = False
    html_context = {}
    print('Warning: Stability of SunPy API page of the documentation requires the ruamel.yaml package to be installed')

# The name of an image file (within the static path) to use as favicon of the
# docs. This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "./logo/favicon.ico"

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = f'{project} v{release}'

# Output file base name for HTML help builder.
htmlhelp_basename = project + 'doc'

# A dictionary of values to pass into the template engine’s context for all pages.
html_context['to_be_indexed'] = ['stable', 'latest']

# -- Options for LaTeX output --------------------------------------------------
# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [('index', project + '.tex', project + ' Documentation', author, 'manual')]

# -- Options for manual page output --------------------------------------------
# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [('index', project.lower(), project + ' Documentation', [author], 1)]

# -- Swap to Napoleon ---------------------------------------------------------
# Remove numpydoc
extensions.remove('numpydoc')
extensions.append('sphinx.ext.napoleon')

# Disable having a separate return type row
napoleon_use_rtype = False
# Disable google style docstrings
napoleon_google_docstring = False

# Extra extensions we want to enable.
extensions += ['sphinx_astropy.ext.edit_on_github', 'sphinx.ext.doctest', 'sphinx.ext.githubpages']

# -- Options for the edit_on_github extension ---------------------------------
# Don't import the module as "version" or it will override the
# "version" configuration parameter
from sunpy import __version__  # noqa isort:skip
edit_on_github_project = "sunpy/sunpy"
if 'dev' not in release:
    edit_on_github_branch = f"{version.split('.')[0]}.{version.split('.')[1]}"
else:
    edit_on_github_branch = "master"
edit_on_github_source_root = ""
edit_on_github_doc_root = "docs"
edit_on_github_skip_regex = '_.*|generated/.*'
github_issues_url = 'https://github.com/sunpy/sunpy/issues/'

# -- Options for the Sphinx gallery -------------------------------------------
if has_sphinx_gallery:
    from sphinx_gallery.sorting import ExplicitOrder
    from sphinx_gallery.sorting import ExampleTitleSortKey
    extensions += ["sphinx_gallery.gen_gallery"]
    path = Path.cwd()
    example_dir = path.parent.joinpath('examples')
    sphinx_gallery_conf = {
        'backreferences_dir': str(path.joinpath('generated', 'modules')),
        'filename_pattern': '^((?!skip_).)*$',
        'examples_dirs': example_dir,
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
        ]),
        'within_subsection_order': ExampleTitleSortKey,
        'gallery_dirs': path.joinpath('generated', 'gallery'),
        'default_thumb_file': path.joinpath('logo', 'sunpy_icon_128x128.png'),
        'abort_on_example_error': False,
        'plot_gallery': True,
        'remove_config_comments': True,
        'doc_module': ('sunpy')
    }

"""
Write the latest changelog into the documentation.
"""
target_file = os.path.abspath("./whatsnew/latest_changelog.txt")
try:
    from sunpy.util.towncrier import generate_changelog_for_docs
    if is_development:
        generate_changelog_for_docs("../", target_file)
except Exception as e:
    print(f"Failed to add changelog to docs with error {e}.")
# Make sure the file exists or else sphinx will complain.
open(target_file, 'a').close()


class Generate(Directive):
    """
    Custom directive to include raw output generated using supplied Python code

    This directive is similar to the ``raw`` directive, except that the raw content is generated
    instead of statically provided.  As with the ``raw`` directive, the required argument specifies
    the format of the output (e.g., ``html`` or ``latex``).

    The optional flag ``html_border`` surrounds HTML output with a black border for visual
    separation.
    """
    required_arguments = 1
    optional_arguments = 0
    option_spec = {'html_border': directives.flag}
    has_content = True

    def run(self):
        # Respect the same disabling options as the ``raw`` directive
        if (not self.state.document.settings.raw_enabled
                or not self.state.document.settings.file_insertion_enabled):
            raise self.warning('"%s" directive disabled.' % self.name)

        attributes = {'format': ' '.join(self.arguments[0].lower().split())}

        source, lineno = self.state_machine.get_source_and_line(self.lineno)

        old_stdout, sys.stdout = sys.stdout, StringIO()
        try:
            # Exceute the Python code and capture its output
            exec('\n'.join(self.content))
            text = sys.stdout.getvalue()

            # Wrap HTML output in a black border if requested
            if attributes['format'] == 'html' and 'html_border' in self.options:
                text = f"<div style='border:1px solid black; padding:3px'>{text}</div>"

            # Append the output in the same fashion as the ``raw`` directive
            raw_node = nodes.raw('', text, **attributes)
            raw_node.source, raw_node.line = source, lineno
            return [raw_node]
        except Exception as e:
            message = f"Unable to execute Python code at {os.path.basename(source)}:{lineno}"
            return [nodes.error(None, nodes.paragraph(text=message)),
                    nodes.paragraph(text=str(e))]
        finally:
            sys.stdout = old_stdout


class MiniGallery(Directive):
    """
    Custom directive to insert a mini-gallery

    The required argument is the qualified name of the object.  The mini-gallery will be the subset
    of gallery examples that make use of that object (from that specific namespace).  There can be
    more than one object named, separated by spaces.
    """
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True

    def run(self):
        # Respect the same disabling options as the ``raw`` directive
        if (not self.state.document.settings.raw_enabled
                or not self.state.document.settings.file_insertion_enabled):
            raise self.warning('"%s" directive disabled.' % self.name)

        obj_list = self.arguments[0].split()

        # Concatenate the backreferences file(s)
        lines = []
        for obj in obj_list:
            path = os.path.join(os.getcwd(), 'generated', 'modules', f'{obj}.examples')
            lines += (FileInput(source_path=path).readlines())[5:]  # slice removes heading

        # Append the end for the gallery
        lines += ['\n',
                  '.. raw:: html\n',
                  '\n',
                  '    <div class="sphx-glr-clear"></div>\n']
        text = ''.join(lines)

        include_lines = statemachine.string2lines(text, convert_whitespace=True)
        self.state_machine.insert_input(include_lines, path)

        return []


def setup(app):
    if not has_sphinx_gallery:
        import warnings
        warnings.warn('The sphinx_gallery extension is not installed, so the '
                      'gallery will not be built. You will probably see '
                      'additional warnings about undefined references due '
                      'to this.')
    if has_yaml:
        app.connect("source-read", rstjinja)

    # Add the custom directives
    app.add_directive('generate', Generate)
    app.add_directive('minigallery', MiniGallery)

    # The theme conf provides a fix for circle ci redirections
    fix_circleci(app)
