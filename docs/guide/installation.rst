************
Installation
************

Requirements
============

sunpy requires Python 3.7 or higher.

Installing Scientific Python and sunpy
======================================

sunpy is part of the wider ecosystem of scientific Python packages for solar physics.
Therefore a working sunpy installation is more about installing the scientific Python ecosystem than sunpy itself.

If you do not currently have a working scientific Python distribution this guide will set you up with the Miniconda, which makes it easy to install and manage your scientific Python packages.

To install the Miniconda Python distribution follow the instructions at
`here <https://docs.conda.io/en/latest/miniconda.html>`__.
Although Miniconda makes it simple to switch between Python versions, we recommend that new users install the latest Python 3.x version of Miniconda.

The reason we choose Miniconda over Anaconda, is mainly due to the size as Anaconda comes with a full install of packages you probably do not need and this way you have more direct control over what has been installed into your Python virtual environment.
Furthermore, you bypass the need for the conda resolver to sort out your root environment which should make conda faster to use.

Installing sunpy using Miniconda
--------------------------------

To install sunpy launch a system command prompt or the 'Anaconda Prompt' (under Windows).
First configure conda for to add the `conda-forge channel <https://conda-forge.org/>`__::

    conda config --add channels conda-forge
    conda config --set channel_priority strict

and now to install sunpy within the default conda virtual environment::

    $ conda install sunpy

This will install sunpy and every package it needs to function.

.. note::
    We strongly recommend using a `Python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__ or a `conda virtual environment. <https://towardsdatascience.com/getting-started-with-python-environments-using-conda-32e9f2779307>`__

Updating sunpy
--------------

You can update to the latest version by running::

    conda update sunpy

Installing sunpy on top of an existing scientific Python environment
--------------------------------------------------------------------

This section assumes you already have everything setup, whether that be conda or a Python virtual environment.
These commands are to be executed within these environments.

Conda
^^^^^

If you want to install sunpy within a pre-existing conda environment, you will want to activate the virtual environment and run::

    $ conda activate <name e.g., base>
    $ conda install sunpy

This assumes you have the conda-forge channel added (as above).

Pip
^^^

This is for installing sunpy within a scientific Python distribution or environment, where ``pip`` has been used to install packages.

To acquire a fully working sunpy installation, simply run::

    pip install sunpy[all]

.. note::
    If this does not work, it could be due to a missing C compiler (e.g., ``gcc`` or ``clang``) that is required to build sunpy at install.
    Getting the compiler either from your system package manager, XCode or Anaconda should address this.

If you have a reason to want a more minimal installation, you can install sunpy with no optional dependencies, however this means a lot of submodules will not import::

    pip install sunpy

It is possible to select which "extra" dependencies you want to install, if you know you only need certain submodules::

    pip install sunpy[map,timeseries]

The available options are: ``[asdf]``, ``[dask]``, ``[database]``, ``[image]``, ``[jpeg2000]``, ``[map]``, ``[net]``, ``[timeseries]``, ``[visualization]``.

If you want to develop sunpy we would strongly recommend reading the `Newcomers' Guide <https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html>`__.

.. note::
    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.

    Do **not** install sunpy or other third-party packages using ``sudo``.

    This error implies you have an incorrectly configured virtual environment or it is not activated.

    If you really do not want to use any virtual environment, you can always do ``pip install --user sunpy``.

Testing sunpy
=============

sunpy provides a method to run the basic test suite that will check that the install has worked correctly.

To run the basic test suite and ensure that your sunpy install is working correctly, use the :func:`sunpy.self_test`::

    import sunpy
    sunpy.self_test()

You will see something like the following in your terminal::

    Starting sunpy self test...
    Checking for packages needed to run sunpy:
    All required and optional sunpy dependencies are installed.
    Starting the sunpy test suite:
    ...

    The tests will run and will report any fails.  You can report these through the `sunpy issue tracker <https://github.com/sunpy/sunpy/issues>`__ and we will strive to help.

It is possible to run this command in a situation where not all packages are installed. If this is the case, you will see the following when you run the test suite::

    Starting sunpy self test...
    Checking for packages needed to run sunpy:
    The following packages are not installed for the sunpy[database] requirement:
    * sqlalchemy
    ...
    You do not have all the required dependencies installed to run the sunpy test suite.
    If you want to run the sunpy tests install the 'tests' extra with `pip install sunpy[all,tests]`

This does not mean sunpy is broken, but you will need to install the extra packages to ensure a "complete" installation of sunpy and run the entire test suite.
It is quite likely that you will run into not having the tests dependencies installed.
