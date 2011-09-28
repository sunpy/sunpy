===============
Latest Unstable
===============
You can install the latest "nightly build" version of SunPy from `GitHub <http://github.com//sunpy/sunpy/>`__.
GitHub provides a convienient download link. It may be easier in the long run to grab SunPy 
using `Git <http://git-scm.com/download>`__, a light and easy version control system: ::

    git clone git://git@github.com/sunpy/sunpy.git

In either case you can now use `paver <http://paver.github.com/>`__ to create a link to the SunPy 
directory you just downloaded in a location that is accessible to Python: ::

    sudo paver develop
    
This will make it make possible to import and use SunPy from anywhere on your system.
Alternatively, if you would prefer you can also skip the paver step and simply
include SunPy locally by starting an interactive Python session from the root SunPy directory.