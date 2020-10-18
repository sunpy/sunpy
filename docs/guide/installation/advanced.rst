.. _advanced-install:

Advanced Installation Instructions
**********************************

This document provides details on how to install and manage your own scientific Python + sunpy installation.
If you have never installed/used a scientific Python environment, we recommend that you follow the :ref:`Miniconda installation instructions <main-install>`.

Installing sunpy on top of an existing scientific Python environment
====================================================================

These instructions assume you have a scientific Python distribution with access to the ``pip`` command installed.

.. warning::
    Users of the Miniconda python distribution should follow the instructions for :ref:`main-install`.
    Mixing packages from conda and pip is not recommended.

Prerequisites
-------------

You will need a C compiler (e.g., ``gcc`` or ``clang``) to be installed as we have a C library within sunpy that is built at install time.
We do provide compiled binaries for sunpy for Linux and Mac OS X (we do not compile our C extension on Windows), so you might not need a C compiler.
If you use Miniconda, you can get these compilers from there.
On Linux, using the package manager for your distribution will usually be the easiest route, while on MacOS X you will need the XCode command line tools.

Using ``pip``
-------------

.. note::
    We strongly recommend using a `Python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__ so any Python packages you install will be kept separate from your Operating System Python.
    You can use Conda as a virtual environment manager instead of one of the more common Python virtual environment managers.

sunpy consists of many submodules that each have their own requirements.
You do not need to fulfil all the requirements if you only intend on using parts of sunpy.
It is however *strongly* recommended to have all the dependencies installed (with the potential exception of ``glymur``).

There are multiple options depending on how many optional dependencies you want to install:

To install sunpy with ``pip`` including optional dependencies (recommended), simply run::

    pip install sunpy[all]

To install sunpy with no optional dependencies::

    pip install sunpy

To install sunpy with net-based dependencies::

    pip install sunpy[net]

To install sunpy with database dependencies::

    pip install sunpy[database]

Other available options are: ``[image]``, ``[jpeg2000]``, ``[asdf]``, ``[tests]`` and ``[docs]``.

The entire list of options are encompassed by a ``[dev]`` option, so you can do::

    pip install sunpy[dev]

to install all the packages needed to run and develop sunpy.
However if you want to develop sunpy we would strongly recommend reading our `newcommers guide <https://docs.sunpy.org/en/latest/dev_guide/newcomers.html>`__.

.. note::
    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.
    In this case you may consider using the ``--user`` option to install the package into your home directory.
    You can read more about how to do this in the `pip documentation <https://pip.pypa.io/en/stable/user_guide/#user-installs>`__.

    Do **not** install sunpy or other third-party packages using ``sudo``.

sunpy's Requirements
********************

sunpy has the following strict requirements:

- `Python <https://www.python.org/>`__ 3.6.x or later.

- `NumPy <https://www.numpy.org/>`__  1.15.0 or later.

- `SciPy <https://www.scipy.org/>`__ 1.0.0 or later.

- `Astropy <https://www.astropy.org/>`__ 3.2.0 or later.

- `pandas <https://pandas.pydata.org/>`__ 0.23.0 or later.

- `parfive <https://pypi.org/project/parfive/>`__ 1.0 or later.

- `matplotlib <https://matplotlib.org/>`__ 2.2.2 or later.

These packages that will be installed as dependencies by default and are the ones required to import the core datatypes `~sunpy.map`, `~sunpy.timeseries` and `~sunpy.spectra`.

sunpy also depends on other packages for optional features.
However, note that these only need to be installed if those you have requested them when you pip install them.
They are installed by default if you use the conda-forge sunpy package.

The following optional packages are:

- `sqlalchemy <https://www.sqlalchemy.org>`__: For the `~sunpy.database` package.

- `scikit-image <https://scikit-image.org/>`__: For `~sunpy.image`.

- `glymur <https://glymur.readthedocs.io/en/latest/>`_ 0.5.9 or later: To enable reading of JPEG2000 files.
  Glymur requires the installation of the `OpenJPEG C library <https://www.openjpeg.org/>`__.

- `beautifulsoup4 <https://www.crummy.com/software/BeautifulSoup/>`_: For `~sunpy.net`.

- `drms <https://pypi.org/project/drms/>`__: For `~sunpy.net`.

- `python-dateutil <https://dateutil.readthedocs.io/en/stable/>`__: For `~sunpy.net`.

- `zeep <https://python-zeep.readthedocs.io/en/master/>`__: For `~sunpy.net`.

- `tqdm <https://github.com/tqdm/tqdm>`__: For `~sunpy.net`.

- `asdf <https://pypi.org/project/asdf/>`__: For `~sunpy.io.special`.

To run the tests:

- `tox <https://tox.readthedocs.io/>`__.

- `hypothesis <https://github.com/HypothesisWorks/hypothesis-python>`__.

- `pytest-astropy <https://github.com/astropy/pytest-astropy>`__.

- `pytest-cov <https://github.com/pytest-dev/pytest-cov>`__.

- `pytest-mock <https://github.com/pytest-dev/pytest-mock>`__.
