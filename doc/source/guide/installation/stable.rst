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
`pip <http://pypi.python.org/pypi/pip>`__, e.g.: ::

    easy_install sunpy
    
or: ::

    pip install sunpy
    
Both of these will download the latest stable version of SunPy and install
it to a standard location for Python libraries on your system.

Using setup.py
--------------
To perform a "traditional" setup.py install, begin by downloading the most 
recent tarball from `the SunPy Downloads Page <http://www.sunpy.org/download/>`__
and extracting it to a suitable location: ::

    tar xzvpf sunpy-<version>.tar.gz
    
Enter the directory you extracted sunpy to, and run: ::

    python setup.py build
    sudo python setup.py install
    
If everything went well, SunPy should now be installed on your system.
