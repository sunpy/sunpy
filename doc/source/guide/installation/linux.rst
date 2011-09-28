=====
Linux
=====

Overview
--------
**Required**

For its basic functioning, SunPy requires several libraries:

* `NumPy <http://numpy.scipy.org/>`__
* `Matplotlib <http://matplotlib.sourceforge.net/>`__
* `PyFITS <http://www.stsci.edu/resources/software_hardware/pyfits>`_

**Optional**

In addition to the required libraries listed above, there are a couple other
optional dependencies which are only needed for certain features.

* `PyQt4 <http://www.riverbankcomputing.co.uk/software/pyqt/download>`__ (SunPy Plotman)
* `Suds <https://fedorahosted.org/suds/>`__ (VSO/HEK support)

In general, it is very easy to install Python and the prereqs needed for SunPy
on Linux systems. Almost all versions of Linux ship with a recent enough version
of Python, so it is unlikely that you will need to do install Python yourself.

For the prereqs, many of them can be installed using 
`pip <http://www.pip-installer.org/en/latest/index.html>`__ or 
`easy_install <http://pypi.python.org/pypi/setuptools>`__, or you can simply 
install them from your distributions software repository using.

Ubuntu
------
On Ubuntu, all of the prereqs can be install using :command:`apt-get`: ::

    sudo apt-get install python-numpy python-matplotlib python-pyfits python-qt4 python-scipy python-suds python-pip git-core ipython

The above command will install the recommended set of libraries and tools 
including both the required and optional dependencies, and also IPython and Git.



