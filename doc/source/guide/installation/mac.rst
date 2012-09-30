========
Mac OS X
========

Because of the wide variety of methods for setting up a Python environment on
Mac OS X, users have a number of options with how to install SunPy and its
dependencies.

Installation using EPD
======================

The simplest method and the recommended one is to install the `Enthought Python
Distribution (EPD) <http://www.enthought.com/products/epd_free.php/>`_, which
includes both Python and a number of the dependencies needed by SunPy. 

Enthought Python Distribution
-----------------------------
The Enthought Python Distribution (EPD) is metapackage which handles the
installation not only of Python, but also a number of other libraries used
by SunPy including NumPy, SciPy, and Matplotlib.

Install EPD
^^^^^^^^^^^
There are a number versions of EPD available, many of which are not free. For
our purposes, however, either of the free versions below will work: ::

 `EPD Free <http://www.enthought.com/products/epd_free.php/>`_
 `EPD Academic <<http://www.enthought.com/products/edudownload.php>>`_

Download whichever of the above versions of EPD is appropriate, and follow the 
`instructions provided on the EPD website <http://www.enthought.com/products/epdgetstart.php?platform=mac>`_ 
to install EPD on your system.

Install pip
^^^^^^^^^^^
Most Python distributions ship with a tool called `easy_install <http://pypi.python.org/pypi/setuptools>`_ 
which assists with installing Python packages.

Although `easy_install`_ is capable of installing most of the dependencies 
needed for SunPy itself, a more powerful tool called `pip <http://pypi.python.org/pypi/pip>`_
provides a more flexible installation (including support for uninstalling, 
upgrading, and installing from remote sources such as GitHub) and should be 
used instead.

Use `easy_install`_ to install `pip`: ::

 sudo easy_install pip

Install remaining dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Now that pip is installed, we can install the remaining libraries needed by
SunPy.

If you installed the EPD Academic version, run: ::

 sudo pip install --upgrade suds
 
If you installed EPD Free version, run: ::

 sudo pip install --upgrade suds
 sudo pip install --upgrade pyfits
 sudo pip install --upgrade pandas
 
Additionally, if you plan to help with SunPy development, some additional 
dependencies are required: ::

 pip install --upgrade pytest pylint paver tox sphinx numpydoc

All done! You are now ready to :doc:`install SunPy itself <index>`.

Other installation methods
==========================

For users who wish to have more control over the installation of Python, several
alternative installation methods are provided, including instructions for
`Macports <http://www.macports.org/>`_ and `Homebrew <http://mxcl.github.com/homebrew/>`_.
These later methods use the `Python.org <http://python.org/>`_ version of 
Python as the basis installation.

Basic
-----
Installation using DMG installers... (TO BE WRITTEN)

Macports
--------

If you use MacPorts to manage your Python distribution, you can simply type
e.g.::

    sudo port install py27-sunpy

to install SunPy. This will automatically install all of the required
dependencies. You may also install SunPy for Python 2.6 (``py26-sunpy``). If
a new version of SunPy is released, you can update your installation using
e.g::

    sudo port selfupdate
    sudo port upgrade py27-sunpy

For more information on using MacPorts to manage your Python installation,
see the following `page <http://astrofrog.github.com/macports-python/>`_.

Homebrew
--------

`Homebrew <http://mxcl.github.com/homebrew/>`_ is a tool for helping to automate
the installation of a number of useful tools and libraries on Mac OS X. It is
similar to Macports, but attempts to improve on some of the Pitfalls of 
Macports.

Note that if you have already installed either fink or Macports on your system,
it is recommended that you uninstall them before using Homebrew.

Python
^^^^^^
To begin, download and install the `latest 2.x version of Python <http://python.org/download/>`_
using the Mac OS X installer available from `Python.org <http://python.org/>`_.

Homebrew
^^^^^^^^
Next, install and update homebrew: ::

 /usr/bin/ruby -e "$(/usr/bin/curl -fsSL https://raw.github.com/mxcl/homebrew/master/Library/Contributions/install_homebrew.rb)"
 brew doctor

Using homebrew, install Qt and some of the other dependencies needed for 
compilation later on by pip: ::

 brew -v install gfortran pkgconfig git openjpeg readline

Pip
^^^
Most Python distributions ship with a tool called `easy_install <http://pypi.python.org/pypi/setuptools>`_ 
which assists with installing Python packages.

Although `easy_install`_ is capable of installing most of
the dependencies needed for SunPy itself, a more powerful tool called 
`pip <http://pypi.python.org/pypi/pip>`__ provides a more flexible installation 
(including support for uninstalling, upgrading, and installing from remote 
sources such as GitHub) and should be used instead.

To begin, use `easy_install`_ to install `pip`: ::

 sudo easy_install pip

Use to install the remaining SunPy dependencies: ::

 sudo pip install --upgrade distribute
 sudo pip install --upgrade ipython
 sudo pip install --upgrade numpy
 sudo pip install --upgrade scipy
 sudo pip install --upgrade pyfits
 sudo pip install --upgrade suds
 sudo pip install --upgrade pandas
 sudo pip install --upgrade matplotlib
 
Additionally, if you plan to help with SunPy development, some additional 
dependencies are required: ::

 pip install --upgrade pytest pylint paver tox sphinx numpydoc

All done!

Trouble-shooting
================




