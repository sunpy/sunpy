************
Installation
************

If you are new to Python and scientific Python then see the :ref:`main-install` section to get setup with the whole environment.
If you already have a working Python / Scientific Python environment then you can skip to the :ref:`advanced-install` section.

Requirements
============

- python >= 3.6
- astropy >= 3.2

sunpy v1.X or greater is compatible with Python 3.6+.
If you need a Python 2 compatible version, you have to install SunPy v0.9.X.

.. _main-install:

Installing Scientific Python and sunpy
======================================

sunpy is part of the wider ecosystem of scientific Python packages for solar physics.
Therefore a working sunpy installation is more about installing the scientific Python ecosystem than sunpy itself.

If you do not currently have a working scientific Python distribution this guide will set you up with the Anaconda
scientific Python distribution, which makes it easy to install and manage your scientific Python packages.
Other options are listed later in the :ref:`advanced-install` section.

To install the Anaconda Python distribution follow the instructions at
`here <https://docs.anaconda.com/anaconda/install/>`__.
Although Anaconda makes it simple to switch between Python versions, we recommend that new users install
the latest Python 3.x version of Anaconda as SunPy v1.X only supports Python 3.6+.

Installing sunpy using Anaconda
-------------------------------

To install SunPy launch a system command prompt or the 'Anaconda Command Prompt' (under Windows).
First configure conda for sunpy downloads::

    conda config --add channels conda-forge

to install sunpy::

    conda install sunpy

You now have a working sunpy installation.
You can check your sunpy install by following the instructions in :ref:`testing-sunpy`.

Updating sunpy
--------------

You can update to the latest version by running::

    conda update sunpy


Advanced sunpy Installation
===========================

If you do not wish to use Anaconda to install Scientific Python or you already have a scientific Python
installation there are other options for installing SunPy:

.. toctree::
  :maxdepth: 2

  advanced.rst
