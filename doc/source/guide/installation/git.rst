====================
Latest Unstable (Git)
====================
To install the latest version of SunPy from GitHub, use Git: ::

    git clone git://git@github.com/sunpy/sunpy.git
    
Next, use `paver <http://paver.github.com/>`__ to create a link to the SunPy 
directory you just downloaded in a location that is accessible to Python: ::

    sudo paver develop
    
Done! You should now be able to import SunPy from anywhere on your system.
Alternatively, if you would prefer you can also skip the paver step and simply
include sunpy locally by starting Python from the root SunPy directory.