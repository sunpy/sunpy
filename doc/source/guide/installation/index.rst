============
Installation
============

Below are instructions for installation SunPy and its prerequisites on 
different platforms.

Installing Prerequisites
------------------------
For instructions on installing the prereqs required by SunPy, choose your OS 
from the list below.

.. toctree::
   :maxdepth: 1
   
   mac
   linux
   win
   
Installing SunPy
----------------
First, you will need to install Python and various prerequisites as described in
the sections above.  You will also need to download and install `Git <http://git-scm.com/download>`__.

There are several different ways to install SunPy, including a quick method and a
developer method.  Depending on your setup, you may need to preface each of the
``pip ...`` commands with ``sudo pip ...``.

Quick installation
^^^^^^^^^^^^^^^^^^
The easiest way to install SunPy is to use pip: ::

    pip install git+https://github.com/sunpy/sunpy.git 
   
This will download and install the latest version of SunPy. To upgrade SunPy at
a later date, you can run: ::

    pip install --upgrade git+https://github.com/sunpy/sunpy.git
    
That's it!
    
Developer installation
^^^^^^^^^^^^^^^^^^^^^^
If you are considering contributing to the development of SunPy, you will likely
want to keep the SunPy code tree in a convenient location.

Open a terminal and cd to a directory where you wish to download SunPy, and 
run: ::

    git clone https://github.com/sunpy/sunpy.git
    
This will download the latest version of SunPy. Finally, cd into the SunPy code
tree and run: ::

    pip install -e .
    
This will make it make possible to import and use SunPy regardless of your
working directory. Now, whenever you pull changes to the Git branch using
``git pull``, your Python installation will automatically see the latest version
of SunPy.
   
Testing your installation
-------------------------

Now you can test your installation. Open a new Python shell and type these 
commands: ::

>>> import sunpy
>>> sunpy.make_map(sunpy.AIA_171_IMAGE).show()

If all goes well you should see an AIA 171 image on your screen.
