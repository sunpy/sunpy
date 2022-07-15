.. _installing:

************
Installation
************

Requirements
============

sunpy requires Python 3.8 or higher.

Installing Scientific Python and ``sunpy``
==========================================

``sunpy`` is part of the wider ecosystem of scientific Python packages for solar physics.
Therefore a working ``sunpy`` installation is more about installing the scientific Python ecosystem than ``sunpy`` itself.

If you do not currently have a working scientific Python distribution this guide will set you up with the Miniforge, which makes it easy to install and manage your scientific Python packages.

`To install the Miniforge Python distribution follow these instructions <https://github.com/conda-forge/miniforge#install>`__.

We do not recommend that you use Anaconda, there are known conflicts with the "defaults" channel and the "conda-forge" channel which has the ``sunpy`` package.

Installing ``sunpy`` using Miniforge
------------------------------------

To install ``sunpy``, launch a system command prompt or terminal.

.. note::

    We strongly recommend using a `Python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__ or a `conda virtual environment. <https://towardsdatascience.com/getting-started-with-python-environments-using-conda-32e9f2779307>`__

Now to install ``sunpy`` within the your virtual environment::

    $ conda install sunpy

This will install ``sunpy`` and each of its dependencies.

Updating ``sunpy``
------------------

You can update to the latest version by running::

    conda update sunpy

Installing ``sunpy`` on top of an existing scientific Python environment
------------------------------------------------------------------------

This section assumes you already have everything setup, whether that be conda or a Python virtual environment.
These commands are to be executed within these environments.

conda
^^^^^

If you want to install ``sunpy`` within a pre-existing conda environment, you will want to activate the virtual environment and run::

    $ conda activate <name e.g., sunpy>
    $ conda install sunpy

This assumes you have the conda-forge channel added (either via Miniforge or manually added).

pip
^^^

This is for installing ``sunpy`` within a scientific Python distribution or environment, where ``pip`` has been used to install packages.

To acquire a fully working ``sunpy`` installation::

    pip install "sunpy[all]"

.. note::
    If this does not work, it could be due to a missing C compiler (e.g., ``gcc`` or ``clang``) that is required to build sunpy at install.
    Getting the compiler either from your system package manager, XCode or Miniforge should address this.

If you have a reason to want a more minimal installation, you can install sunpy with no optional dependencies, however this means a lot of submodules will not import::

    pip install "sunpy"

It is possible to select which "extra" dependencies you want to install, if you know you only need certain submodules::

    pip install "sunpy[map,timeseries]"

The available options are: ``[asdf]``, ``[dask]``, ``[database]``, ``[image]``, ``[jpeg2000]``, ``[map]``, ``[net]``, ``[timeseries]``, ``[visualization]``.

If you want to develop ``sunpy`` we would strongly recommend reading the `Newcomers' Guide <https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html>`__.

.. note::
    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.

    Do **not** install ``sunpy`` or other third-party packages using ``sudo``.

    This error implies you have an incorrectly configured virtual environment or it is not activated.
