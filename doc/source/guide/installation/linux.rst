=====
Linux
=====

Pre-requisites
--------------
To begin, install the pre-requisites for SunPy. 

Ubuntu
^^^^^^
On Ubuntu, all of the prereqs can be install using :command:`apt-get`: ::

    sudo apt-get install python-numpy python-matplotlib python-pyfits python-qt4 python-scipy python-suds git-core ipython 

The above command will install the recommended set of libraries and tools including both the required and optional dependencies.

Installing SunPy
----------------
There are several different ways to install the latest stable version of SunPy.

Using easy_install or pip
^^^^^^^^^^^^^^^^^^^^^^^^^
The easist way to install SunPy is to use either 
`easy_install <http://peak.telecommunity.com/DevCenter/EasyInstall>`__ or `pip <http://pypi.python.org/pypi/pip>`__, e.g.: ::

    sudo easy_install sunpy
    
or: ::

    sudo pip install sunpy

Using setup.py
^^^^^^^^^^^^^^
To perform a "traditional" setup.py install, begin by downloading the most 
recent tarball from `the SunPy Downloads Page <http://www.sunpy.org/download/>`__
and extracting it to a suitable location: ::

    tar xzvpf sunpy-<version>.tar.gz
    
Enter the directory you extracted sunpy to, and run: ::

    python setup.py build
    sudo python setup.py install
    
If everything went well, SunPy should now be installed on your system.


Installing SunPy Using Git
--------------------------
To install the latest version of SunPy from GitHub, use Git: ::

    git clone git://git@github.com/sunpy/sunpy.git
    
Next, use `paver <http://paver.github.com/>`__ to create a link to the SunPy 
directory you just downloaded in a location that is accessible to Python:

    sudo paver develop
    
Done! You should now be able to import SunPy from anywhere on your system.
Alternatively, if you would prefer you can also skip the paver step and simply
include sunpy locally by starting Python from the root SunPy directory.


Testing your installation
-------------------------

Now you can test your installation. Open a new Python shell by typing 
:command:`python` in ``Command Prompt``, and type these commands: ::

>>> import sunpy
>>> sunpy.Map(sunpy.AIA_171_IMAGE).show()