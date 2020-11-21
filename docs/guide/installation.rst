************
Installation
************

Requirements
============

These are the minimum versions of packages needed to install sunpy:

- python 3.7
- astropy 4.0.3
- scipy 1.2
- parfive 1.1

Installing Scientific Python and sunpy
======================================

sunpy is part of the wider ecosystem of scientific Python packages for solar physics.
Therefore a working sunpy installation is more about installing the scientific Python ecosystem than sunpy itself.

If you do not currently have a working scientific Python distribution this guide will set you up with the Miniconda, which makes it easy to install and manage your scientific Python packages.

To install the Miniconda Python distribution follow the instructions at
`here <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`__.
Although Miniconda makes it simple to switch between Python versions, we recommend that new users install
the latest Python 3.x version of Miniconda as SunPy only supports Python 3.7+.

The reason we choose Miniconda over Anaconda, is mainly due to the size as Anaconda comes with a full install of packages you probably do not need and this way you have more direct control over what has been installed into your Python virtual environment.
Furthermore, you bypass the need for the conda resolver to sort out your root environment which should make "conda" faster to use.

We strongly recommend using a `Python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__ so any Python packages you install will be kept separate from your Operating System Python.
You can use Miniconda as a virtual environment manager instead of one of the more common Python virtual environment managers and that is what the following seciton is about.

Installing sunpy using Miniconda
--------------------------------

To install SunPy launch a system command prompt or the 'Anaconda Prompt' (under Windows).
First configure conda for to add the `conda-forge channel <https://conda-forge.org/>__`::

    conda config --add channels conda-forge
    conda config --set channel_priority strict

and now to install sunpy within a new virtual environment::

    $ conda create -n sunpy python sunpy
    $ conda activate sunpy

This will install sunpy and every package it needs to function.

Updating sunpy
--------------

You can update to the latest version by running::

    conda update sunpy

Installing sunpy on top of an existing scientific Python environment
--------------------------------------------------------------------

Conda
^^^^^

If you want to install sunpy within a pre-existing conda environment, you can run::

    $ conda install sunpy

This assumes you have the conda-forge channel added (as above).

Pip
^^^

This is for installing sunpy within a scientific Python distribution, where ``pip`` has been used to install  packages.

We do provide compiled binaries for sunpy for Linux and Mac OS X (we do not compile our C extension on Windows), so you might not need a C compiler.
If there is no compiled binary, you will need a C compiler (e.g., ``gcc`` or ``clang``) to be installed as we have a C library within sunpy that is built at install time.
If you use Miniconda, you can get these compilers from there.
On Linux, using the package manager for your distribution will usually be the easiest route, while on MacOS X you will need the XCode command line tools.

To acquire a fully working sunpy installation, simply run::

    pip install sunpy[all]

If you have a reason to bypass this, you can sunpy with no optional dependencies::

    pip install sunpy

To install sunpy with map- and timeseries-based dependencies::

    pip install sunpy[map,timeseries]

The available options are: ``[asdf]``, ``[dask]``,``[database]``,``[instr]``,``[image]``,``[jpeg2000]``,``[map]``,``[net]``,``[timeseries]``,``[visualization]``.

If you want to develop sunpy we would strongly recommend reading the `Newcomers' Guide <https://docs.sunpy.org/en/latest/dev_guide/newcomers.html>`__.

.. note::
    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.

    Do **not** install sunpy or other third-party packages using ``sudo``.

    This error implies you either an incorrectly configured virtual environment or it was not activate.

    If you really do not want to use any virtual environment, you can always do ``pip install --user sunpy``.
