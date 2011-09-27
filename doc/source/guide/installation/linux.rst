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

