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
>>> sunpy.Map(sunpy.AIA_171_IMAGE).show()

Mac
---
Installation instructions for Mac.

Windows
-------
There are many ways to get SunPy up and running on Windows, and we describe the recommended method and an alternate method below.  Lines in text boxes should be typed into ``Command Prompt``, which can be started from the Start menu or in Start->Run by typing :command:`cmd`.

Recommended method
^^^^^^^^^^^^^^^^^^

**1. Install Python(x,y)**

Download and install `Python(x,y) <https://code.google.com/p/pythonxy/wiki/Downloads>`_.  Python(x,y) is a distribution that include not only Python, but also a large variety of Python modules and development tools.  Please note that this installer is rather large (>500 MB) and thus may take a while to download.

**2. Install additional modules**

Install PyFITS, which is used to read FITS files.  (On some configurations, there currently may be a bug where you will not be able to read compressed image data from FITS files.  We are working to resolve this bug.) ::

    easy_install pyfits

Install Paver, which is used to link symbolically to the SunPy code. ::

    easy_install paver

**3. Install Git**

Download and install `Git <https://code.google.com/p/msysgit/downloads/list?can=3>`_ (choose the first file listed).  Git is used to retrieve the SunPy code.

**4. Download and install SunPy**

The following will download SunPy to ``C:\sunpy``.  If you wish to download SunPy elsewhere, modify these and later commands accordingly. ::

    cd C:\
    "%ProgramFiles%\Git\bin\git" clone git://github.com/sunpy/sunpy.git

If you get the error ``The system cannot find the path specified``, try using these commands: ::

    cd C:\
    "%ProgramFiles(x86)%\Git\bin\git" clone git://github.com/sunpy/sunpy.git

Now that SunPy is downloaded, you can create a symbolic link for Python to find the SunPy code. ::

    cd C:\sunpy\
    paver develop

In the future, to update SunPy to the latest version: ::

    cd C:\sunpy\
    "%ProgramFiles%\Git\bin\git" pull

As before, if you get the error ``The system cannot find the path specified``, try using these commands: ::

    cd C:\sunpy\
    "%ProgramFiles(x86)%\Git\bin\git" pull


Alternate method
^^^^^^^^^^^^^^^^

Please use this method only if you are experienced with computers and cannot use the recommended method.  Possible reasons include having very little disk space or needing to have the most up-to-date versions of modules.

**1a. Install Python**

Download and install `Python 2.7 <http://www.python.org/ftp/python/2.7.2/python-2.7.2.msi>`_ (32-bit).  Note that even if you are on a 64-bit system, you are installing a 32-bit version of Python to be able to use precompiled binaries.

You should update the ``PATH`` environment variable so that Python executables and associated scripts can be found:

    1. Go to ``Start``-> ``Control Panel`` -> ``System`` -> ``Advanced system settings`` -> ``Environment variables``
    2. Find the ``PATH`` environment variable, under either user or system variables, and append ``C:\Python27`` and ``C:\Python27\Scripts``, separated by semicolons.
    

**1b. Install packaged modules**

Download and install `NumPy <http://sourceforge.net/projects/numpy/files/NumPy/1.6.1/numpy-1.6.1-win32-superpack-python2.7.exe/download>`_.

Download and install `SciPy <http://sourceforge.net/projects/scipy/files/scipy/0.9.0/scipy-0.9.0-win32-superpack-python2.7.exe/download>`_.

Download and install `matplotlib <http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.0.1/matplotlib-1.0.1.win32-py2.7.exe/download>`_.

Download and install `setuptools 
<http://pypi.python.org/packages/2.7/s/setuptools/setuptools-0.6c11.win32-py2.7.exe>`_.


**2-4. The remaining steps**

You have now performed the required elements of step 1 of the recommended method.  Now perform steps 2-4 of that method to complete your installation.


Test your installation
^^^^^^^^^^^^^^^^^^^^^^

Now you can test your installation. Open a new Python shell by typing :command:`python` in ``Command Prompt``, and type these commands: ::

>>> import sunpy
>>> sunpy.Map(sunpy.AIA_171_IMAGE).show()
