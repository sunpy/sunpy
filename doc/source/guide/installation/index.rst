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
There are several different ways to install SunPy. For most users the first
method will be sufficient. If you would like to become more involved in the
project, then an alternative method using Git is recommended.

Quick installation
^^^^^^^^^^^^^^^^^^
Once you have successfully installed Python and the prereqs as described in the
sections above you are ready to install SunPy itself.

The easiest way to install SunPy is to use pip: ::

    sudo pip install git+https://github.com/sunpy/sunpy.git 
   
This will download and install the latest version of SunPy. To upgrade the
installation when a new version comes out you can run the same command with
the '--upgrade' switch: ::

    sudo pip install --upgrade git+https://github.com/sunpy/sunpy.git
    
That's it!
    
Developer installation
^^^^^^^^^^^^^^^^^^^^^^
If you are considering contributing to the development of SunPy, a more flexible
approach using `Git <http://git-scm.com/>`__ and `Paver <http://paver.github.com/>`__
is suggested.

The first step is to `download and install Git <http://git-scm.com/download>`__.

Next, we will use pip to install Paver -- a Python tool similar to `GNU Make <http://www.gnu.org/software/make/>`__: ::

    sudo pip install paver
    
Now, open a terminal and cd to a directory where you wish to download SunPy, and 
run: ::

    git clone https://github.com/sunpy/sunpy.git
    
This will download the latest version of SunPy. Finally, run the following
paver command: ::

    sudo paver develop
    
This will make it make possible to import and use SunPy from anywhere on your 
system. Now, whenever you pull changes to the Git branch using ``git pull``, your SunPy
installation will be automatically up to date.
   
Testing your installation
-------------------------

Now you can test your installation. Open a new Python shell and type these 
commands: ::

>>> import sunpy
>>> sunpy.make_map(sunpy.AIA_171_IMAGE).show()

If all goes well you should see an AIA 171 image on your screen.
