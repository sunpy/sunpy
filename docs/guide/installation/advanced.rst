.. _advanced-install:

Advanced Installation Instructions
**********************************

This document provides details on things you need to know to install and manage your own scientific Python + SunPy installation.
If you have never installed or used scientific Python we recommend that you follow the :ref:`Anaconda installation instructions <main-install>`.

Installing SunPy on top of an existing Scientific Python Environment
====================================================================

These instructions assume you have a scientific Python distribution with access to the `pip` command installed.

Prerequisites
-------------

You will need a compiler suite and the development headers for Python and numpy in order to build SunPy.
On Linux, using the package manager for your distribution will usually be the easiest route, while on MacOS X you will need the XCode command line tools.

The `instructions for building Numpy from source <https://docs.scipy.org/doc/numpy/user/install.html>`_ are also a good resource for setting up your environment to build Python packages.

Using `pip`
-----------

.. warning::
    Users of the Anaconda python distribution should follow the instructions for :ref:`anaconda_install`.

SunPy consists of many submodules that each have their own requirements.
You do not need to fulfil all the requirements if you only intend on using parts of SunPy.
It is however *strongly* recommended to have all the dependencies installed (with the potential exception of ``glymur``).

There are multiple options depending on how many optional dependencies you want to install:

To install SunPy with ``pip`` including optional dependencies (recommended), simply run::

    pip install sunpy[all]

To install SunPy with no optional dependencies::

    pip install sunpy

To install SunPy with net-based dependencies::

    pip install sunpy[net]

To install SunPy with database dependencies::

    pip install sunpy[database]

Other available options are: ``[image]``, ``[jpeg2000]``, ``[asdf]``, ``[tests]`` and ``[docs]``.

The entire list of options are compassed by a ``[dev]`` option, so you can do::

    pip install sunpy[dev]

to install all the packages needed to run and develop SunPy.

.. note::
    You will need a C compiler (e.g., ``gcc`` or ``clang``) to be installed.
    If you use anaconda, you can get these packages from there instead.

.. note::
    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.
    In this case you may consider using the ``--user`` option to install the package into your home directory.
    You can read more about how to do this in the `pip documentation <https://pip.pypa.io/en/stable/user_guide/#user-installs>`__.

    Alternatively, if you intend to do development on other software that uses SunPy, such as an affiliated package, consider installing SunPy into a `virtualenv <https://docs.python-guide.org/dev/virtualenvs/>`__.

    Do **not** install SunPy or other third-party packages using ``sudo`` unless you are fully aware of the risks.

.. _testing-sunpy:

Testing SunPy
*************

.. warning::
    The tests will fail if you do not install all the optional dependencies.

The easiest way to test your installed version of SunPy is running correctly is to use the :func:`sunpy.self_test()` function::

    import sunpy
    sunpy.self_test()

which will run many of the SunPy tests.

The tests should run and print out any failures, which you can report at the `SunPy issue tracker <https://github.com/sunpy/sunpy/issues>`__.

SunPy's Requirements
********************

SunPy has the following strict requirements:

- `Python <https://www.python.org/>`__ 3.6.x or later.

- `NumPy <https://www.numpy.org/>`__  1.15.0 or later.

- `SciPy <https://www.scipy.org/>`__ 1.0.0 or later.

- `Astropy <https://www.astropy.org/>`__ 3.2.0 or later.

- `pandas <https://pandas.pydata.org/>`__ 0.23.0 or later.

- `parfive <https://pypi.org/project/parfive/>`__ 1.0 or later.

- `matplotlib <https://matplotlib.org/>`__ 2.2.2 or later.

SunPy also depends on other packages for optional features.
However, note that these only need to be installed if those particular features are needed.
SunPy will import even if these dependencies are not installed.

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

The packages that will be installed as dependencies by default and are the ones
required to import the core datatypes `~sunpy.map`, `~sunpy.timeseries` and
`~sunpy.spectra`. These are the strict requirements and the following optional
packages:
