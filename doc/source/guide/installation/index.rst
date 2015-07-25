============
Installation
============

SunPy is a Python package for solar physics, it relies on and enables the use
of the wider ecosystem of scientific Python packages for solar physics.
Therefore getting a working SunPy installation is more about installing the
whole scientific Python ecosystem than SunPy itself.

If you are new to Python and scientific Python then continue to follow this
guide to get setup with the whole environment. If you already have a working
Python / Scientific Python environment then you can skip to the
:ref:`advanced-install` section.

.. _main-install:

Installing Scientific Python and SunPy
--------------------------------------

If you do not currently have a working scientific Python distribution this
guide will set you up with the Anaconda Python distribution. Alternative options
exist for different operating systems and can be found later in the
:ref:`advanced-install` section.

Anaconda is a free distribution of Python and a large amount of common
scientific packages, it is very powerful and easy to use. Installing Anaconda
provides (almost) all the packages you need to use SunPy.

To install the Anaconda Python distribution follow the instructions
`here <http://docs.continuum.io/anaconda/install.html>`_. You will need to
select the correct download for your platform and follow the install procedure.

Installing SunPy on top of Anaconda
###################################

To install SunPy launch a system command prompt or the
'Anaconda Command Prompt' (under Windows). First configure conda
for sunpy downloads::

    conda config --add channels sunpy --add channels astropy

to install SunPy::

    conda install sunpy

You now have a working SunPy installation. You can check your SunPy install
by following the instructions in :ref:`testing-sunpy`.

.. note::

    Currently Glymur / JPEG2000 support is not tested under Anaconda on any
    platforms. However, an external installation of openJPEG should work with
    the glymur installation in conda.

Updating SunPy to a New Version
###############################

When a new version of SunPy is released you can update to the latest version
by running::

    conda update sunpy


Advanced SunPy Installation
---------------------------

If you do not wish to use Anaconda to install Scientific Python or you already
have a scientific Python installation there are other options for installing
SunPy.

.. toctree::
  :maxdepth: 2

  advanced.rst
