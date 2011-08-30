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
Installation instructions for Windows.  Lines in text boxes should be typed into ``Command Prompt``, which can be started from the Start menu or in Start->Run by typing ``cmd``.

Recommended method
^^^^^^^^^^^^^^^^^^

**1. Install Python(x,y)**

Download and install `Python(x,y) <https://code.google.com/p/pythonxy/wiki/Downloads>`_.  Python(x,y) is a distribution that include not only Python, but also a large variety of Python modules and development tools.  Please note that this installer is rather large (>500 MB) and thus may take a while to download.  (An alternative distribution is the `Enthought Python Distribution (free version) <http://www.enthought.com/products/epd_free.php>`_.)

**2. Install additional modules**

Install PyFITS, which is used to read FITS files.  It is highly likely you will see a warning message that an optional extension module ``pyfits.pyfitsComp`` failed to build.  If you do not require the ability to read compressed image data from FITS files, then you can ignore this message.  Otherwise, see the section below.

    cd C:\Python26\Scripts\
    easy_install pyfits

Install Paver, which is used to link symbolically to the SunPy code

    cd C:\Python26\Scripts\
    easy_install paver

**3. Install Git**

Download and install `Git <https://code.google.com/p/msysgit/downloads/list?can=3>`_.  Git is used to retrieve the SunPy code.

**4. Download and install SunPy

The following assumes that you will download SunPy to ``C:\sunpy``.  If you wish to download SunPy elsewhere, modify these commands accordingly.

    cd C:\
    %ProgramFiles%\Git\bin\git clone git://github.com/sunpy/sunpy.git
    cd C:\sunpy\
    paver develop

In the future, to update SunPy to the latest version:

    cd C:\sunpy\
    %ProgramFiles%\Git\bin\git pull


Alternate method
^^^^^^^^^^^^^^^^

Please use this method only if you are experienced with computers and for some reason wish to install only the bare minimum of packages to run SunPy.

**1a. Install Python**

Download and install the latest version of Python 2.7 from the `official page <http://www.python.org/getit/>`_.

Next, in order to enable Python to find libraries installed on your computer
you need to update the ``PATH`` environmental variable on your machine:

    1. Click ``Start``-> ``Control Panel`` -> ``System`` -> ``Advanced system settings`` -> ``Environment variables``
    2. Find the ``PATH`` environmental variable under either user or system variables and the filepath to your Python installation (e.g. "C:\Python27").
    

**1b. Install packaged modules**

Download and install `NumPy <http://sourceforge.net/projects/numpy/files/NumPy/1.6.1/numpy-1.6.1-win32-superpack-python2.7.exe/download>`_.

Download and install `SciPy <http://sourceforge.net/projects/scipy/files/scipy/0.9.0/scipy-0.9.0-win32-superpack-python2.7.exe/download>`_.

Download and install `matplotlib <http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.0.1/matplotlib-1.0.1.win32-py2.7.exe/download>`_.

Download and install `setuptools 
<http://pypi.python.org/packages/2.7/s/setuptools/setuptools-0.6c11.win32-py2.7.exe>`_.


**2–4. The remaining steps**

You have now performed the equivalent of step 1 of the recommended method.  Perform steps 2–4 of that method to complete your installation.  The only difference is that you should type ``C:\Python27\`` in commands rather than ``C:\Python26\``


Test your installation
^^^^^^^^^^^^^^^^^^^^^^

To test it all out, open a new Python shell and try typing: ::

>>> import sunpy
>>> sunpy.Map(sunpy.AIA_171_IMAGE).plot()


"pyfits.pyfitsComp failed to build"?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This error (and any others like it) are likely due to Python being unable to find a compatible C++ compiler in your path.  Resolving this error may not be for the faint-of-heart.

Our suggested compiler is the `MinGW <http://www.mingw.org/>`_ compiler.  If you have followed the "recommended method" of installation, then you already have such a compiler (MinGW), and in fact it is already in your path.  However, Python is not configured to use MinGW by default.

To configure Python to use MinGW by default, create the file ``C:\Python26\lib\distutils\distutils.cfg`` containing these lines:

    [build]
    compiler=mingw32
    [build_ext]
    compiler=mingw32

