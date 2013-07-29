=======
Windows
=======

Overview
--------

There are many ways to get the necessary environment up and running for SunPy on Windows, 
and we describe the recommended method and an alternate method below.  

Python(x,y)
^^^^^^^^^^^

**1. Install Python(x,y)**

Download and install `Python(x,y) <https://code.google.com/p/pythonxy/wiki/Downloads>`_.
Python(x,y) is a distribution that include not only Python, but also a large 
variety of Python modules and development tools.  Please note that this 
installer is rather large (~600 MB) and thus may take a while to download.
The modules in Python(x,y) include:

* NumPy
* SciPy
* Matplotlib
* PyFITS
* pandas
* distribute
* pip
* Qt
* Sphinx
* pylint
* ipython

**2. Set up the Python build environment**

Create a file in ``C:\Python27\lib\distutils\`` called ``distutils.cfg`` that
contains the following lines: ::

    [build]
    compiler=mingw32
    [build_ext]
    compiler=mingw32

**3. Install Git**

Download and install `Git <https://code.google.com/p/msysgit/downloads/list?can=3&q=Full+installer+for+official+Git+for+Windows>`_.
Git is used to retrieve the SunPy code.

You are now ready to continue with
:doc:`the rest of the installation <index>`.

Alternate method
^^^^^^^^^^^^^^^^

If Python(x,y) is unavailable, or if it is critical for you to have the most
up-to-date version of modules, then this alternate method eschews installers
as much as possible and emphasizes building from source.  This method is not
for the faint of heart!  Please use this method only if you are experienced
with computers and Python.

**1. Install Python**

Download and install `Python 2.7 <http://www.python.org/ftp/python/2.7.5/python-2.7.5.msi>`_ 
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

Unfortunately, the compilations of SciPy and Matplotlib will likely fail due to
missing libraries.  Until there is a workable solution, you should download the
latest installers: `SciPy <http://sourceforge.net/projects/scipy/files/scipy/0.12.0/scipy-0.12.0-win32-superpack-python2.7.exe/download>`_
and `Matplotlib <http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.2.1/matplotlib-1.2.1.win32-py2.7.exe/download>`__.
You are now ready to continue with :doc:`the rest of the installation <index>`.

