=============
Latest Stable
=============
Overview
--------

To install the latest stable version of SunPy, choose one of the methods listed
below. In general, it is recommended that you use `pip` for easier upgrading 
and removal.

Using easy_install or pip
-------------------------
The easist way to install SunPy is to use either 
`easy_install <http://peak.telecommunity.com/DevCenter/EasyInstall>`__ or 
`pip <http://pypi.python.org/pypi/pip>`__, (Mac users please refer to this `page <http://sunpy.org/doc/guide/installation/mac.html#installing-other-packages>`__
for more instructions) e.g.: ::

    easy_install sunpy
    
or: ::

    pip install sunpy
    
Both of these will download the latest stable version of SunPy and install
it to a standard location for Python libraries on your system.

Using setup.py
--------------
To perform a "traditional" setup.py install, begin by downloading the most 
recent version from `the SunPy Downloads Page <http://www.sunpy.org/download/>`__
Enter the directory you extracted Sunpy to, and run: ::

    python setup.py build
    sudo python setup.py install
    
(Mac users please refer to this `page <http://sunpy.org/doc/guide/installation/mac.html#installing-other-packages>`__
for more instructions). If everything went well, SunPy should now be installed on your system.
