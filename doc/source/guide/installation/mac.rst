========
Mac OS X
========

Overview
--------

Because of the wide variety of methods for setting up a Python environment on
Mac OS X, users have a number of options with how to install SunPy and its
dependencies. We have a number of solutions listed below. Choose the solution which best
suits you.

Install Anaconda
----------------

To install the Anaconda Python distribution follow the instructions `here <http://docs.continuum.io/anaconda/install.html>`_.

.. warning::
    
    You must have the XCode developer package installed to be able to compile SunPy.
    see :ref:`xcode`.

Once complete follow the install instructions for the :ref:`anaconda_install`.

Other installation methods
--------------------------

For users who wish to have more control over the installation of Python, several
alternative installation methods are provided below, including instructions for
`Macports <http://www.macports.org/>`_ and `Homebrew <http://mxcl.github.com/homebrew/>`_.
The following instructions are not recommended for beginners. OS X comes pre-loaded with
Python but each versions of OS X (Mountain Lion, Snow Leopard, etc.) ships with a
different version of Python. In the instructions below, we therefore recommend that you
install your own version of Python. You will then have two versions of Python living on
your system at the same time. It can be confusing to make sure that when you install
packages they are installed for the correct Python so beware.

Python 2.7.3+
^^^^^^^^^^^^^
Download and install the latest version of 
`32 bit Python 2.7.3 <http://www.python.org/download/releases/2.7.3/>` 
using their DMG installer. Next, choose your package installer of choice (either
Macports or Homebrew) and follow the instructions below. If you do not have either
go to their respective websites and install one of the other as needed.

.. _xcode:

XCode tools / Compiler
^^^^^^^^^^^^^^^^^^^^^^
If you are using MacOS X, you will need to the XCode command line
tools.  One way to get them is to install `XCode
<https://developer.apple.com/xcode/>`__. If you are using OS X 10.7
(Lion) or later, you must also explicitly install the command line
tools. You can do this by opening the XCode application, going to
**Preferences**, then **Downloads**, and then under **Components**,
click on the Install button to the right of **Command Line Tools**.
Alternatively, on 10.7 (Lion) or later, you do not need to install
XCode, you can download just the command line tools from
https://developer.apple.com/downloads/index.action (requires an Apple
developer account).
 
Macports
^^^^^^^^

Macports is a useful To install everything that you need using Macports. First install Macports. You'll
need Xcode and the command line tools. First update the package index: ::

    sudo port selfupdate

To install the basic Python packages, run: ::

    sudo port install py27-matplotlib py27-numpy py27-scipy 

To install Qt (this will take a while!), run: ::

    sudo port install qt4-mac

This will install python 2.7 at the same time. To make the macports version of Python
the default type the following into a terminal: ::

    sudo port select --set python python27

For more information on using MacPorts to manage your Python installation,
see the following `page <http://astrofrog.github.com/macports-python/>`_.

Homebrew
^^^^^^^^

`Homebrew <http://mxcl.github.com/homebrew/>`_ is a tool for helping to automate
the installation of a number of useful tools and libraries on Mac OS X. It is
similar to Macports, but attempts to improve on some of the pitfalls of 
Macports.

Note that if you have already installed either fink or Macports on your system,
it is recommended that you uninstall them before using Homebrew.

Next, install and update homebrew: ::

 /usr/bin/ruby -e "$(/usr/bin/curl -fsSL https://raw.github.com/mxcl/homebrew/master/Library/Contributions/install_homebrew.rb)"
 brew doctor

Using homebrew, install Qt and some of the other dependencies needed for 
compilation later on by pip: ::

 brew -v install gfortran pkgconfig git openjpeg readline pyqt

Now on to the next steps.

Git
^^^
Head over and `download <http://git-scm.com/downloads>` and install git. 

Pip
^^^
Most Python distributions ship with a tool called 
`easy_install <http://pypi.python.org/pypi/setuptools>`_ 
which assists with installing Python packages.

Although `easy_install`_ is capable of installing most of
the dependencies needed for SunPy itself, a more powerful tool called 
`pip <http://pypi.python.org/pypi/pip>`__ which provides a more flexible installation 
(including support for uninstalling, upgrading, and installing from remote 
sources such as GitHub) and should be used instead. 

Use `easy_install`_ to install `pip`: ::

    sudo easy_install pip

You are now ready to install scipy, numpy, and matplotlib.

Scientific Libraries
^^^^^^^^^^^^^^^^^^^^
If pip installed properly, then you can install NumPy simply with: ::

    pip install numpy
    
Now under Lion, install the stable version of SciPy (0.10) by running: ::

    pip install scipy

Mountain Lion users will need to install the development version of SciPy (0.11) 
by executing the following line:

pip install -e git+https://github.com/scipy/scipy#egg=scipy-dev

Now on to matplotlib

On Lion, install matplotlib like any other package: ::

    pip install matplotlib

Mountain Lion users will have to use the development version as of this writing: ::

    pip install git+https://github.com/matplotlib/matplotlib.git#egg=matplotlib-dev

Done! You are now ready to :doc:`install SunPy itself <index>`.

Installation using Virtual Box
------------------------------

As a last resort you can always install Linux in a Virtual Machine.
`Virtual box <https://www.virtualbox.org/>` is a free virtual environment that allows 
you to run linux or other operating systems concurrently with OS X. Since it is very easy
to install and maintain SunPy under linux, the idea here is to have your own linux 
environment on your mac and use it for SunPy. We recommend you use Ubuntu linux in your 
virtual box. You can download an install disk for Ubuntu at 
`their website <http://www.ubuntu.com/download/help/install-desktop-latest>`. Directions
on how to install your first virtual machine on Virtual Box are also 
`available <https://www.virtualbox.org/manual/ch01.html#gui-createvm>`. After your have
Ubuntu installed just follow the :doc:`instructions for linux <linux>`! Simple.

