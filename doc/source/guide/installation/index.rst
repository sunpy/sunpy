============
Installation
============

Below are instructions for installation SunPy and its prerequisites on 
different platforms.

Installing Scientific Python
----------------------------
For instructions on setting up the scientific Python environment which is required by SunPy, 
choose your OS from the list below. When you are done come back here and follow
the instructions below to install SunPy and its requirements.

.. toctree::
   :maxdepth: 1
   
   mac
   linux
   win

Installing Python Modules
-------------------------
You should now have the following  on your system; Python, Numpy, Scipy, 
Matplotlib, git, pip and Qt. You can now install some final remaining dependencies using pip.
Depending on which system you are using some of these may already be installed 
but it does not hurt to upgrade: ::

 pip install --upgrade distribute
 pip install --upgrade pyfits
 pip install --upgrade suds
 pip install --upgrade pandas
 pip install --upgrade beautifulsoup4

We also recommend you use `ipython <http://ipython.org/>` (a great python enviromnent). 
This can also be installed using pip: ::

 pip install --upgrade ipython
 
Additionally, if you plan to help with SunPy development, some additional 
dependencies are required: ::

 pip install --upgrade pytest pylint paver tox sphinx numpydoc

All done with the SunPy Python prerequisites. You are now ready to install SunPy itself.

Installing SunPy
----------------
There are a number of ways of installing SunPy. We currently recommend you use the latest
(development) version of SunPy as we are developing quickly. To do this you'll use git 
and pip to grab the latest code from github. If you would rather have a more stable
version of SunPy then check out the section below on grabbing the stable SunPy.

Grabbing the latest SunPy
=========================

There are two ways to install the latest version SunPy, including a quick method 
and a developer method.  Depending on your setup, you may need to preface each of the
``pip ...`` commands with ``sudo pip ...``.

Quick installation
^^^^^^^^^^^^^^^^^^
The easiest way to install SunPy is to use pip (combined with git): ::

    pip install git+https://github.com/sunpy/sunpy.git 
   
This will download and install the latest version of SunPy. To upgrade SunPy at
a later date, you can run: ::

    pip install --upgrade --no-deps git+https://github.com/sunpy/sunpy.git

Please make sure to include the ``--no-deps`` option because otherwise pip may
try to upgrade dependencies such as SciPy and Matplotlib that are difficult to
build from source and the likely errors will abort the upgrade. That's it!
    
Developer installation
^^^^^^^^^^^^^^^^^^^^^^
If you are considering contributing to the development of SunPy, you will likely
want to keep the SunPy code tree in a convenient location.

Open a terminal and cd to a directory where you wish to download SunPy, and 
run: ::

    git clone https://github.com/sunpy/sunpy.git
    
This will download the latest version of SunPy. Finally, from inside the new SunPy 
directory run: ::

    pip install -e .
    
This will make it make possible to import and use SunPy regardless of your
working directory. To upgrade SunPy at a later date, go into the SunPy directory 
and run: ::

    git pull upstream master
    
That's it!
   
Grabbing the stable SunPy
=========================
If you'd like to use a more stable version of SunPy and only upgrade at milestones
then head over to `github <https://github.com/sunpy/sunpy/tags>` and download the latest
stable release. After extracting the download you can install SunPy by using pip like so
: ::

    pip install ./downloads/SunPyPackage-1.0.4.tar.gz
    
This will add it to your existing Python modules. To upgrade to the next version, just
head back to the github page and grab the next version and repeat the process above.

Testing your installation
-------------------------

Now you can test your installation. Open a new Python shell and type these 
commands: ::

>>> import sunpy
>>> sunpy.make_map(sunpy.AIA_171_IMAGE).peek()

If all goes well you should see an AIA 171 image on your screen.
