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

Installing sunpy using Miniconda
--------------------------------

To install SunPy launch a system command prompt or the 'Anaconda Prompt' (under Windows).
First configure conda for to add the `conda-forge channel <https://conda-forge.org/>__`::

    conda config --add channels conda-forge
    conda config --set channel_priority strict

and now to install sunpy::

    conda install sunpy

This will install sunpy and every package it needs to function.

Updating sunpy
--------------

You can update to the latest version by running::

    conda update sunpy

Installing sunpy on top of an existing scientific Python environment
--------------------------------------------------------------------

These instructions assume you have a scientific Python distribution with access to the ``pip`` command installed.

Prerequisites
^^^^^^^^^^^^^

We do provide compiled binaries for sunpy for Linux and Mac OS X (we do not compile our C extension on Windows), so you might not need a C compiler.
If there is no compiled binary, you will need a C compiler (e.g., ``gcc`` or ``clang``) to be installed as we have a C library within sunpy that is built at install time.
If you use Miniconda, you can get these compilers from there.
On Linux, using the package manager for your distribution will usually be the easiest route, while on MacOS X you will need the XCode command line tools.

Using ``pip``
^^^^^^^^^^^^^

.. note::
    We strongly recommend using a `Python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__ so any Python packages you install will be kept separate from your Operating System Python.
    You can use Conda as a virtual environment manager instead of one of the more common Python virtual environment managers.

sunpy consists of many sub-modules that each have their own requirements.
You do not need to fulfil all the requirements if you only intend on using parts of sunpy.
It is however it **strongly recommended** to have all the dependencies installed.
This ensure that your version of sunpy is fully operational.

To install sunpy with ``pip`` including all optional dependencies (**recommended**), simply run::

    pip install sunpy[all]

To install sunpy with no optional dependencies::

    pip install sunpy

To install sunpy with map- and timeseries-based dependencies::

    pip install sunpy[map,timeseries]

The available options are: ``[asdf]``, ``[dask]``,``[database]``,``[instr]``,``[image]``,``[jpeg2000]``,``[map]``,``[net]``,``[timeseries]``,``[visualization]``.

If you want to develop sunpy we would strongly recommend reading the `Newcomers' Guide <https://docs.sunpy.org/en/latest/dev_guide/newcomers.html>`__.

.. note::
    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.
    In this case you may consider using the ``--user`` option to install the package into your home directory.
    You can read more about how to do this in the `pip documentation <https://pip.pypa.io/en/stable/user_guide/#user-installs>`__.
    Furthermore, if this does occur, it is most likely that your virtual environment is incorrect configured.

    Do **not** install sunpy or other third-party packages using ``sudo``.
