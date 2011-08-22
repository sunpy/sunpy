------------
Installation
------------
Pre-requisites
--------------
SunPy stands on the shoulders of giants:

* `NumPy <http://numpy.scipy.org/>`_
* `Matplotlib <http://matplotlib.sourceforge.net/>`_
* `PyFITS <http://www.stsci.edu/resources/software_hardware/pyfits>`_

Linux
-----
Installation instructions for Linux.

Ubuntu
^^^^^^
To begin, install the pre-requisites for SunPy using :command:`apt-get`: ::

    sudo apt-get install python-numpy python-matplotlib python-pyfits python-scipy git-core ipython

The ``ipython`` package in the above list installs the `IPython enhanced console 
<http://ipython.scipy.org/moin/>`_ and is optional but recommended.

Next, use Git to download a copy of the latest version of SunPy: ::

    git clone git://git@github.com/sunpy/sunpy.git

Done! To see if everything went okay, start a Python session and try importing
SunPy:

>>> import sunpy
>>> sunpy.Map(sunpy.AIA_171_IMAGE).plot()

Mac
---
Installation instructions for Mac.

Windows
-------
Installation instructions for Windows.

Recommended method
^^^^^^^^^^^^^^^^^^

**1. Install Python(x,y)**

**2. Install additional modules**

Finally, SunPy uses the PyFITS library to read 
`FITS <http://en.wikipedia.org/wiki/FITS>`_ files. PyFITS does
not offer a Windows binrary, so we will need to install it some other way.
Fortunately, there is a great tool for Python called ``easy_install`` which 
makes installing many libraries very trivial. EasyInstall is part of the
Setuptools package.



Install paver

**4. Install Git**

TODO

To begin, open a command-line shell and run the following command:

    git clone git://github.com/sunpy/sunpy.git

This will download the latest version of SunPy. You then need to copy the 
folder sunpy (from inside the root sunpy directory) to ``C:\Python27\Lib`` so 
that Python can find it.


Alternate method
^^^^^^^^^^^^^^^^

**1. Install Python**

To begin, grab the latest version of the Python 2.x for Windows from the
`Python downloads page <http://www.python.org/getit/>`_.  Run the installer
and follow the instructions on screen.


Next, in order to enable Python to find libraries installed on your computer
you need to update the ``PATH`` environmental variable on your machine:

    1. Click ``Start``-> ``Control Panel`` -> ``System`` -> ``Advanced system settings`` -> ``Environment variables``
    2. Find the ``PATH`` environmental variable under either user or system variables and the filepath to your Python installation (e.g. "C:\Python27").
    

**2. Install additional modules**

Next we will install `NumPy <http://numpy.scipy.org/>`_, `SciPy 
<http://www.scipy.org/>`_, `Matplotlib <http://matplotlib.sourceforge.net/>`_, 
and `IPython <http://ipython.scipy.org/moin/>`_. Install each package by running
the corresponding downloaded executable.  They should  automatically find the 
location of the Python installation.

    1. To begin, download and install the latest `NumPy Windows binary <http://sourceforge.net/projects/numpy/files/NumPy/1.6.0b2/numpy-1.6.0b2-win32-superpack-python2.7.exe/download>`_.
    2. Next, download and install the latest version of the `SciPy Windows binary <http://sourceforge.net/projects/scipy/files/scipy/0.9.0/scipy-0.9.0-win32-superpack-python2.7.exe/download>`_
    3. Download and install the latest version of `Matplotlib for Windows <http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.0.1/matplotlib-1.0.1.win32-py2.7.exe/download>`_.
    4. Finally, as an option (but recommended step), download and install `IPython for Windows <http://ipython.scipy.org/dist/0.10.1/ipython-0.10.1.win32-setup.exe>`_.
    

To install `Setuptools 
<http://pypi.python.org/pypi/setuptools>`_, simply download and run the `latest
windows installer 
<http://pypi.python.org/packages/2.7/s/setuptools/setuptools-0.6c11.win32-py2.7.exe>`_.

To install PyFITS, open a console and call ``easy_install.exe pyfits`` from 
inside the ``scripts`` directory of your Python installation: ::

    cd C:\Python27\Scripts
    easy_install.exe pyfits

**3. Continue with recommended approach**


Test installation
^^^^^^^^^^^^^^^^^

To test it all out, open a new Python shell and try typing: ::

>>> import sunpy
>>> sunpy.Map(sunpy.AIA_171_IMAGE).plot()


C++ compiler
^^^^^^^^^^^^

MinGW
distutils.cfg
