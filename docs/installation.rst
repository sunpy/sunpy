.. _installing:

************
Installation
************

The SunPy project maintains a range of libraries that leverage the wider ecosystem of scientific Python packages for solar physics.

Installing a SunPy library
==========================

If you do not currently have a working scientific Python distribution this guide will set you up with `Miniforge <https://conda-forge.org/#about>`__, which makes it easy to install and manage Python packages.

`To install the Miniforge Python distribution follow these instructions <https://github.com/conda-forge/miniforge#install>`__.

We do not recommend that you use Anaconda, there are known conflicts with the "defaults" channel and the "conda-forge" channel that contains all of the SunPy libraries.

Installing a package using Miniforge
------------------------------------

To install ``sunpy`` (or another package such as ``ablog``), launch a terminal (under a UNIX-like system) or Miniforge Prompt (under Windows).

.. note::

    We strongly recommend using a `Python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__ or a `conda virtual environment. <https://towardsdatascience.com/getting-started-with-python-environments-using-conda-32e9f2779307>`__

You will want to activate your virtual environment then run::

    $ conda install sunpy

This will install ``sunpy`` and each of its dependencies.

Updating a package
------------------

You can update to the latest version of any package by running::

    conda update <package_name>

Installing on top of an existing scientific Python environment
--------------------------------------------------------------

This section assumes you already have everything setup, whether that be conda or a Python virtual environment.
These commands are to be executed within these environments.

conda
^^^^^

If you want to install ``sunpy`` within a pre-existing conda environment, you will want to activate the virtual environment and run::

    $ conda install sunpy

This assumes that your pre-existing conda environment is already using the "conda-forge" channel.
If this is not the case, please install Miniforge (using the instructions above).

pip
^^^

This is for installing ``sunpy`` within a scientific Python distribution or environment, where ``pip`` has been used to install packages.

To acquire a fully working ``sunpy`` installation::

    pip install "sunpy[all]"

.. note::

    We strive to provide binary wheels for all of our packages.
    If you are using a Python distribution or operating system that is missing a binary wheel.
    ``pip`` will try to compile the package from source and this is likely to fail without a C compiler (e.g., ``gcc`` or ``clang``).
    Getting the compiler either from your system package manager or XCode should address this.

If you have a reason to want a more minimal installation, you can install sunpy with no optional dependencies, however this means a lot of submodules will not import::

    pip install "sunpy"

It is possible to select which "extra" dependencies you want to install, if you know you only need certain submodules::

    pip install "sunpy[map,timeseries]"

The available options are: ``[asdf]``, ``[dask]``, ``[database]``, ``[image]``, ``[jpeg2000]``, ``[map]``, ``[net]``, ``[timeseries]``, ``[visualization]``.

If you want to develop ``sunpy`` we would strongly recommend reading the `Newcomers' Guide <https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html>`__.

.. note::
    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.
    Do **not** install ``sunpy`` or other Python package using ``sudo``.
    This error implies you have an incorrectly configured virtual environment or it is not activated.
