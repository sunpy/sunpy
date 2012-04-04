=====
Linux
=====

Overview
--------
**Required**

For its basic functioning, SunPy requires several libraries:

* `NumPy <http://numpy.scipy.org/>`__
* `SciPy <http://www.scipy.org/>`__
* `Matplotlib <http://matplotlib.sourceforge.net/>`__ (1.0+)
* `PyFITS <http://www.stsci.edu/resources/software_hardware/pyfits>`_
* `Suds <https://fedorahosted.org/suds/>`__
* `pandas <http://pandas.sourceforge.net/dsintro.html>`_
* `PyQt4 <http://www.riverbankcomputing.co.uk/software/pyqt/download>`__

**Optional**

In addition to the required libraries listed above, if you plan to work with
JPEG 2000 data, you must also install:

* `OpenJPEG <http://www.openjpeg.org/>`__
* `PIL <http://www.pythonware.com/products/pil/>`__

In general, it is very easy to install Python and the prereqs needed for SunPy
on Linux systems. Almost all versions of Linux ship with a recent enough version
of Python, so it is unlikely that you will need to do install Python yourself.

For the prereqs, many of them can be installed using 
`pip <http://www.pip-installer.org/en/latest/index.html>`__ or 
`easy_install <http://pypi.python.org/pypi/setuptools>`__, or you can simply 
install them from your distributions software repository using.

Ubuntu
------
On Ubuntu, most of the pre-reqs are available in the Ubuntu software repos and
can be installed using :command:`apt-get`: ::

    sudo apt-get install python-numpy python-matplotlib python-pyfits python-qt4 python-scipy python-suds python-imaging python-pip openjpeg-tools git-core ipython

The above command will install the recommended set of libraries and tools 
including both the required and optional dependencies, and also IPython and Git.

Next, use pip to install pandas:

    sudo pip install pandas



