.. _advanced-install:

Advanced Installation Instructions
----------------------------------

This document provides details on things you need to know to install and
manage your own scientific Python + SunPy installation.
If you have never installed or used scientific Python we recommend that you
follow the :ref:`Anaconda installation instructions <main-install>`.


Alternative Scientific Python Installation Instructions
=======================================================

There are many alternatives to Anaconda as a way of installing a scientific
Python environment, there are various other platform specific ways to install
scientific Python:

.. toctree::
    :maxdepth: 1

    linux.rst
    mac.rst


Installing SunPy on top of an existing Scientific Python Environment
====================================================================

These instructions assume you have a scientific Python distribution with access
to the `pip` command installed.

Prerequisites
#############

You will need a compiler suite and the development headers for Python and
Numpy in order to build SunPy. On Linux, using the package manager for your
distribution will usually be the easiest route, while on MacOS X you will
need the XCode command line tools.

The `instructions for building Numpy from source
<http://docs.scipy.org/doc/numpy/user/install.html>`_ are also a good
resource for setting up your environment to build Python packages.

.. note:: If you are using MacOS X, you will need to the XCode command line
          tools.  One way to get them is to install `XCode
          <https://developer.apple.com/xcode/>`_. If you are using OS X 10.7
          (Lion) or later, you must also explicitly install the command line
          tools. You can do this by opening the XCode application, going to
          **Preferences**, then **Downloads**, and then under **Components**,
          click on the Install button to the right of **Command Line Tools**.
          Alternatively, on 10.7 (Lion) or later, you do not need to install
          XCode, you can download just the command line tools from
          `the Apple developer site <https://developer.apple.com/downloads/index.action>`_
          (requires an Apple developer account).

SunPy's Requirements
####################

SunPy consists of many submodules that each have their own requirements.
You do not need to fulfil all the requirements if you only intend on using
parts of SunPy. It is however *strongly* recommended to have all the
dependencies installed (with the potential exception of `glymur`).

SunPy has the following strict requirements:

- `Python <http://www.python.org/>`_ 2.7

- `NumPy <http://www.numpy.org/>`_  1.6.0 or later

- `SciPy <http://www.scipy.org/>`_ 0.10.0 or later

- `AstroPy <http://www.astropy.org/>`__ 1.0.0 or later

SunPy also depends on other packages for optional features.
However, note that these only need to be installed if those particular features
are needed. SunPy will import even if these dependencies are not installed.

- `Matplotlib <http://matplotlib.org/>`_ [*Highly Recommended*] 1.3.0 or later:
  For `~sunpy.lightcurve`, `~sunpy.map`, `~sunpy.spectra`, `~sunpy.instr`
  and `~sunpy.visualization`.

- `pandas <http://pandas.pydata.org/>`_ 0.10 or later: For `~sunpy.lightcurve`.

- `sqlalchemy <http://www.sqlalchemy.org>`_: For the `~sunpy.database` package.

- `suds <https://bitbucket.org/jurko/suds>`_: For `~sunpy.net`.

- `beautifulsoup4 <http://www.crummy.com/software/BeautifulSoup/>`_: For `~sunpy.spectra.Callisto` Spectrograms and `~sunpy.net.helio`

- `requests <http://docs.python-requests.org/en/latest/>`_: For the `~sunpy.net.jsoc` submodule.

- `wcsaxes <http://wcsaxes.readthedocs.org/en/latest/>`_: For `sunpy.map` plotting improvements.

- `glymur <https://glymur.readthedocs.org/en/latest/>`_ 0.5.9 or later: To enable reading of JPEG2000 files.
  Glymur requires the installation of the `OpenJPEG C library <http://code.google.com/p/openjpeg/downloads/list>`_.

- `pytest <http://pytest.org/latest/>`_: To run tests.

The packages that will be installed as dependencies by default and are the ones
required to import the core datatypes `~sunpy.map`, `~sunpy.lightcurve` and
`~sunpy.spectra`. These are the strict requirements and the following optional
packages:

Using `pip`
###########

There are multiple options depending on how many optional dependencies you
want to install:

To install SunPy with `pip` including optional dependencies (recommended), simply run::

    pip install sunpy[all]

To install SunPy with no optional dependencies::

    pip install sunpy

To install SunPy with net-based dependencies (suds and beautifulsoup)::

    pip install sunpy[net]

To install SunPy with database dependencies (sqlalchemy)::

    pip install sunpy[database]

.. warning::
    Users of the Anaconda python distribution should follow the instructions
    for :ref:`anaconda_install`.

.. note::

    You will need a C compiler (e.g. ``gcc`` or ``clang``) to be installed.

.. note::

    If you get a ``PermissionError`` this means that you do not have the
    required administrative access to install new packages to your Python
    installation.  In this case you may consider using the ``--user`` option
    to install the package into your home directory.  You can read more about
    how to do this in the `pip documentation
    <http://www.pip-installer.org/en/1.2.1/other-tools.html#using-pip-with-the-user-scheme>`_.

    Alternatively, if you intend to do development on other software that uses
    SunPy, such as an affiliated package, consider installing SunPy into a
    :ref:`virtualenv<using-virtualenv>`.

    Do **not** install SunPy or other third-party packages using ``sudo``
    unless you are fully aware of the risks.

.. _testing-sunpy:

Testing SunPy
-------------

The easiest way to test your installed version of SunPy is running
correctly is to use the :func:`sunpy.self_test()` function::

    import sunpy
    sunpy.self_test(online=False)

which will run many of the SunPy tests.

The tests should run and print out any failures, which you can report at
the `SunPy issue tracker <http://github.com/sunpy/sunpy/issues>`_.


Installing the Development Version of SunPy
-------------------------------------------

The latest (bleeding-edge) development version of SunPy can be cloned from github
using this command::

   git clone git://github.com/sunpy/sunpy.git

.. note::

   If you wish to participate in the development of SunPy, see
   :ref:`developer-docs`.  This document covers only the basics
   necessary to install SunPy.

Once inside the source directory that has been clone from GitHub you can install
SunPy using::

    python setup.py install

.. note::
    This command will need access to system folders append `--user` to install
    SunPy into your home directory.

Troubleshooting
---------------

If you get an error mentioning that you do not have the correct permissions to
install SunPy into the default ``site-packages`` directory, you should try
installing with::

    pip install sunpy --user

which will install into a default directory in your home directory.

Building documentation
----------------------

.. note::
    Building the documentation is in general not necessary unless you
    are writing new documentation or do not have internet access, because
    the latest (and archive) versions of SunPy's documentation are
    available at `docs.sunpy.org <http://docs.sunpy.org>`_ .

Building the documentation requires the SunPy source code and some additional
packages:

    - `Sphinx <http://sphinx.pocoo.org>`_ (and its dependencies) 1.0 or later

    - `Graphviz <http://www.graphviz.org>`_

.. note::

    Sphinx also requires a reasonably modern LaTeX installation to render
    equations.  Per the `Sphinx documentation
    <http://sphinx-doc.org/builders.html?highlight=latex#sphinx.builders.latex.LaTeXBuilder>`_,
    for the TexLive distribution the following packages are required to be
    installed:

    * latex-recommended
    * latex-extra
    * fonts-recommended

    For other LaTeX distributions your mileage may vary. To build the PDF
    documentation using LaTeX, the ``fonts-extra`` TexLive package or the
    ``inconsolata`` CTAN package are also required.

There are two ways to build the SunPy documentation. The most straightforward
way is to execute the command (from the sunpy source directory)::

    python setup.py build_sphinx -lo

The documentation will be built in the ``doc/build/html`` directory, and can
be read by pointing a web browser to ``doc/build/html/index.html``.

The LaTeX documentation can be generated by using the command::

    python setup.py build_sphinx -b latex

The LaTeX file ``SunPy.tex`` will be created in the ``doc/build/latex``
directory, and can be compiled using ``pdflatex``.

The above method builds the API documentation from the source code.
