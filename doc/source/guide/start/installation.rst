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
Installation instructions for linux.

Ubuntu
^^^^^^
To begin, install the pre-requisites for SunPy using apt-get: ::

    sudo apt-get install python-numpy python-matplotlib python-pyfits python-scipy bzr ipython

The ``ipython`` package in the above list installs the `IPython enhanced console 
<http://ipython.scipy.org/moin/>`_ and is optional but recommended.

Next, use Bazaar to download a copy of the latest version of SunPy: ::

    bzr branch lp:sunpy

Done! To see if everything went okay, start a Python session and try importing
sunpy:

>>> import sunpy
>>> sunpy.Map(‘sunpy/dev/sample-data/AIA20110319_105400_0171.fits’).plot()




