============
Installation
============

SunPy is a Python package for solar physics, it relies on and enables the use
of the wider ecosystem of scientific Python packages. Therefore getting a 
working SunPy installation is more about installing the whole scientific Python
ecosystem than SunPy itself.

If you are new to Python and scientific Python then continue to follow this 
guide to get setup with the whole environment. If you already have a working
Python / Scientific Python environment then you can skip to the <REQUIREMENTS>
section.

Installing Scientific Python
----------------------------

If you do not currently have a working scientific Python distribution this 
guide will set you up with the Anaconda Python distribution. Alternative options
exist for different operating systems and can be found later in this guide.

Anaconda is a free distribution of Python and a large amount of common 
scientific packages, it is very powerful and easy to use. Installing Anaconda 
provides (almost) all the packages you need to use SunPy.

To install the Anaconda Python distribution follow the instructions 
`here <http://docs.continuum.io/anaconda/install.html>`_. You will need to 
select the correct download for your platform and follow the install procedure.

.. note::

    On OS/X you need to install XCode so you can build SunPy's extensions.
    see :ref:`xcode`

Anaconda comes with it's own package manager (`conda`) for installing 
packages provided by Continuum Analytics, it is strongly recommended to use the
conda command for installing and updating packages where possible.
For packages not availible through `conda` (like SunPy) it is possible to 
install them from the source distribution using the Python `pip` command.

Installing SunPy on top of Anaconda
###################################

To install SunPy launch the Anaconda command prompt or a system prompt and 
run the following commands:

To install SunPy's extra dependancies run::

    conda update astropy
    pip install suds

To install run::
 
 	pip install sunpy

You now have a working SunPy installation, you can now skip to 
<TESTING YOUR SUNPY INSTALL>.

If a new version of SunPy is released you can update to the latest version 
by running::

    conda update anaconda
    pip install --upgrade --no-deps sunpy suds

.. note::

    Currently Glymur / JPEG2000 support is not tested under Anaconda on any 
    platforms. If you require JPEG2000 support either use a different install 
    method, or contact the SunPy mailing list.

