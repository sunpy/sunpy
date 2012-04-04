=======
Windows
=======

Overview
--------
**Required**

For its basic functioning, SunPy requires several libraries:

* `NumPy <http://numpy.scipy.org/>`__
* `SciPy <http://www.scipy.org/>`__
* `Matplotlib <http://matplotlib.sourceforge.net/>`__ (1.0+)
* `PyFITS <http://www.stsci.edu/resources/software_hardware/pyfits>`_
* `Suds <https://fedorahosted.org/suds/>`__
* `pandas <http://pandas.pydata.org/>`_
* `PyQt4 <http://www.riverbankcomputing.co.uk/software/pyqt/download>`__

**Optional**

In addition to the required libraries listed above, if you plan to work with
JPEG 2000 data, you must also install:

* `OpenJPEG <http://www.openjpeg.org/>`__
* `PIL <http://www.pythonware.com/products/pil/>`__

There are many ways to get SunPy up and running on Windows, and we describe the 
recommended method and an alternate method below.  Lines in text boxes should 
be typed into ``Command Prompt``, which can be started from the Start menu or 
in Start->Run by typing :command:`cmd`.

Recommended method
^^^^^^^^^^^^^^^^^^

**1. Install Python(x,y)**

Download and install `Python(x,y) <https://code.google.com/p/pythonxy/wiki/Downloads>`_.
Python(x,y) is a distribution that include not only Python, but also a large 
variety of Python modules and development tools.  Please note that this 
installer is rather large (~400 MB) and thus may take a while to download.

**2. Install other required modules**

Download and install `pip <http://pythonxy.googlecode.com/files/pip-1.0.2_py27.exe>`_.  (Note: this installer is built by the Python(x,y) team.)

Download and install `PyFITS <http://pypi.python.org/packages/2.7/p/pyfits/pyfits-3.0.3.win32-py2.7.exe>`__.

**3. Install optional modules**

Install optional modules: ::

    pip install suds
    pip install paver

* Paver is used to link symbolically to the SunPy code

**4. Install Git**

Download and install `Git <https://code.google.com/p/msysgit/downloads/list?can=3>`_ 
(choose the first file listed).  Git is used to retrieve the SunPy code.

**5. Download and install SunPy**

The following will download SunPy to ``C:\sunpy``.  If you wish to download 
SunPy elsewhere, modify these and later commands accordingly. ::

    cd C:\
    "%ProgramFiles%\Git\bin\git" clone git://github.com/sunpy/sunpy.git

If you get the error ``The system cannot find the path specified``, try using 
these commands: ::

    cd C:\
    "%ProgramFiles(x86)%\Git\bin\git" clone git://github.com/sunpy/sunpy.git

Now that SunPy is downloaded, you can create a symbolic link for Python to find 
the SunPy code. ::

    cd C:\sunpy\
    paver develop

In the future, to update SunPy to the latest version: ::

    cd C:\sunpy\
    "%ProgramFiles%\Git\bin\git" pull

As before, if you get the error ``The system cannot find the path specified``, 
try using these commands: ::

    cd C:\sunpy\
    "%ProgramFiles(x86)%\Git\bin\git" pull


Alternate method
^^^^^^^^^^^^^^^^

Please use this method only if you are experienced with computers and cannot 
use the recommended method.  Possible reasons include having very little disk 
space or needing to have the most up-to-date versions of modules.

**1a. Install Python**

Download and install `Python 2.7 <http://www.python.org/ftp/python/2.7.2/python-2.7.2.msi>`_ 
(32-bit).  Note that even if you are on a 64-bit system, you are installing a 
32-bit version of Python to be able to use precompiled binaries.

You should update the ``PATH`` environment variable so that Python executables 
and associated scripts can be found:

    1. Go to ``Start``-> ``Control Panel`` -> ``System`` -> ``Advanced system settings`` -> ``Environment variables``
    2. Find the ``PATH`` environment variable, under either user or system variables, and append ``C:\Python27`` and ``C:\Python27\Scripts``, separated by semicolons.
    

**1b. Install required and optional modules included in Python(x,y)**

Download and install `NumPy <http://sourceforge.net/projects/numpy/files/NumPy/1.6.1/numpy-1.6.1-win32-superpack-python2.7.exe/download>`__.

Download and install `SciPy <http://sourceforge.net/projects/scipy/files/scipy/0.9.0/scipy-0.9.0-win32-superpack-python2.7.exe/download>`__.

Download and install `matplotlib <http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.0.1/matplotlib-1.0.1.win32-py2.7.exe/download>`__.

Download and install `PyQt4 <http://www.riverbankcomputing.co.uk/static/Downloads/PyQt4/PyQt-Py2.7-x86-gpl-4.8.5-1.exe>`__.

Download and install `distribute <http://pythonxy.googlecode.com/files/distribute-0.6.21_py27.exe>`_.  (Note: this installer is built by the Python(x,y) team.)

**2-5. The remaining steps**

You have now performed the required elements of step 1 of the recommended 
method.  Now perform steps 2-5 of that method to complete your installation.

.. _NumPy: http://numpy.scipy.org/
.. _SciPy: http://www.scipy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net/>
.. _PyFITS: http://www.stsci.edu/resources/software_hardware/pyfits>
