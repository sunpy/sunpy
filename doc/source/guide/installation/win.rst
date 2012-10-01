=======
Windows
=======

Overview
--------

There are many ways to get SunPy up and running on Windows, and we describe the 
recommended method and an alternate method below.  Lines in text boxes should 
be typed into ``Git Bash``, which you can find in your Start Menu after you
install Git.  (Alternatively, you can use ``Command Prompt`` if you know where
the Git binaries are located or you have opted to have them in your path.)

Recommended method
^^^^^^^^^^^^^^^^^^

**1. Install Python(x,y)**

Download and install `Python(x,y) <https://code.google.com/p/pythonxy/wiki/Downloads>`_.
Python(x,y) is a distribution that include not only Python, but also a large 
variety of Python modules and development tools.  Please note that this 
installer is rather large (~500 MB) and thus may take a while to download.

**2. Install other required modules**

The following required modules are already included in Python(x,y): NumPy,
SciPy, Matplotlib, PyFITS, pandas, and distribute.  Download and install
`pip <http://code.google.com/p/pythonxy/downloads/list?q=pip>`_.  (Note: this
installer is built by the Python(x,y) team.)

**3. Install Git**

Download and install `Git <https://code.google.com/p/msysgit/downloads/list?can=3&q=Full+installer+for+official+Git+for+Windows>`_.
Git is used to retrieve the SunPy code.

**4. Install recommended/optional modules**

The following recommended/optional modules are already included in Python(x,y):
IPython, PIL, pylint, and PyQt.  You can use pip to install additional modules,
e.g.: ::

    pip install suds
    pip install beautifulsoup4
    pip install pytest

**5. Download and install SunPy**

The simple (non-developer) way to install SunPy is to use the combination of pip
and Git: ::

    pip install git+https://github.com/sunpy/sunpy.git

On the other hand, if you are a developer and plan to contribute to SunPy, it is
often easier to have the SunPy code in a more accessible location that you
manage with Git.  In that case, you can use pip to install your working version
of the SunPy code (here assumed to be located in ``C:\sunpy``): ::

    pip install -e 'C:\sunpy'

**6. Upgrading in the future**

If you used the simple (non-developer) way to install SunPy, you can upgrade
SunPy by using pip and Git again: ::

    pip install --upgrade --no-deps git+https://github.com/sunpy/sunpy.git

Please make sure to include the ``--no-deps`` option because otherwise pip may
try to upgrade dependencies such as SciPy and Matplotlib that are difficult to
build from source and the likely errors will abort the upgrade.

If you used the developer way to install SunPy, you should use Git as
appropriate to upgrade SunPy.

Python(x,y) and its included modules (listed above in steps 2 and 4) can be
upgraded by downloading and installing the latest versions.  Git can also be
upgraded the same way.  Recommended/optional modules not included in Python(x,y)
can generally be upgraded by using pip, e.g.: ::

    pip install --upgrade suds
    pip install --upgrade beautifulsoup4
    pip install --upgrade pytest


Alternate method
^^^^^^^^^^^^^^^^

If Python(x,y) is unavailable, or if it is critical for you to have the most
up-to-date version of modules, then this alternate method eschews installers
as much as possible and emphasizes building from source.  This method is not
for the faint of heart!  Please use this method only if you are experienced
with computers and Python.

**1. Install Python**

Download and install `Python 2.7 <http://www.python.org/ftp/python/2.7.3/python-2.7.3.msi>`_ 
(32-bit).  Note that even if you are on a 64-bit system, you are installing a 
32-bit version of Python to be able to use precompiled binaries where still needed.

You should update the ``PATH`` environment variable so that Python executables 
and associated scripts can be found:

    1. Go to ``Start``-> ``Control Panel`` -> ``System`` -> ``Advanced system settings`` -> ``Environment variables``
    2. Find the ``PATH`` environment variable, under either user or system variables, and append ``C:\Python27`` and ``C:\Python27\Scripts``, separated by semicolons.

**2. Install compilers**

Download and install `MinGW <http://mingw.org/>`_, specifically the C, C++, and
Fortran compilers.  Make sure that the binaries can be found in your path (e.g.,
by adding ``C:\MinGW\bin`` to your path).

**3. Set up the Python build environment**

Create a file in ``C:\Python27\lib\distutils\`` called ``distutils.cfg`` that
contains the following lines: ::

    [build]
    compiler=mingw32
    [build_ext]
    compiler=mingw32

There is currently a bug in the Python 2.7 code, so you will also need to edit
``cygwincompiler.py`` in the same directory.  Remove all five instances of the
character string "-mno-cygwin".

**4. Install pip**

Download `distribute_setup.py <http://python-distribute.org/distribute_setup.py>`_
to be able to install distribute: ::

    python distribute_setup.py

Download `get-pip.py <https://raw.github.com/pypa/pip/master/contrib/get-pip.py>`_
to be able to install pip: ::

    python get-pip.py

**5. Install required modules**

You can use pip to download and build modules from source: ::

    pip install numpy
    pip install scipy
    pip install matplotlib
    pip install pyfits
    pip install pandas

Unfortunately, the compilations of SciPy and Matplotlib will likely fail due to
missing libraries.  Until there is a workable solution, you should download the
latest installers: `SciPy <http://sourceforge.net/projects/scipy/files/scipy/0.11.0/scipy-0.11.0-win32-superpack-python2.7.exe/download>`_
and `Matplotlib <http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.1.1/matplotlib-1.1.1.win32-py2.7.exe/download>`__.

**6. The remaining steps**

You have completed the essential elements of steps 1-2 of the recommended 
method.  Continue with steps 3-5 of that method to complete your installation.

.. _NumPy: http://numpy.scipy.org/
.. _SciPy: http://www.scipy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net/>
.. _PyFITS: http://www.stsci.edu/resources/software_hardware/pyfits>
