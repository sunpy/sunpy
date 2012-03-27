===============
Latest Unstable
===============
You can install the latest "nightly build" version of SunPy from 
`GitHub <http://github.com/sunpy/sunpy/>`__. GitHub provides a convienient 
download link. It may be easier in the long run to grab SunPy using 
`Git <http://git-scm.com/download>`__, a light and easy version control system: ::

    git clone git://git@github.com/sunpy/sunpy.git
    
Installing with Python 2
------------------------

If you are using Python 2, you can use `paver <http://paver.github.com/>`__ to 
create a link to the SunPy directory you just downloaded in a location that is 
accessible to Python: ::

    sudo paver develop
    
This will make it make possible to import and use SunPy from anywhere on your 
system. Alternatively, if you would prefer you can also skip the paver step 
and simply include SunPy locally by starting an interactive Python session 
from the root SunPy directory.

Installing with Python 3
------------------------
SunPy also supports installation for Python 3. Instead of using `paver develop`,
however, you must instead perform a more routine installation using setup.py.
Additionally, it is recommended that you first install the Python 
`distribute <http://pypi.python.org/pypi/distribute>`__ package before 
installing SunPy.