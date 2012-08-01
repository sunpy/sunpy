=====
Linux
=====

Overview
--------

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

    sudo apt-get update && sudo apt-get install python-numpy  \
    python-matplotlib python-pyfits python-qt4 python-scipy python-suds \
    python-imaging python-pandas python-pip openjpeg-tools git-core ipython

The above command will install the recommended set of libraries and tools 
including both the required and optional dependencies, and also IPython and Git.
    
Arch Linux
----------
The install the latest version of SunPy on Arch Linux, run: ::

    yaourt sunpy

(or some similar command using the AUR installer of your choice...)