============
Installation
============

Below are instructions for installing SunPy and its prerequisites on different 
platforms.

SunPy's Dependancies
--------------------

SunPy consists of many submodules that each have their own requirements. You do
not need to fulfil all the requirements if you only intend on using parts of 
SunPy. It is of course recommended that you install the requirements for all 
of SunPy.

The Full set of requirements is:

* `NumPy <http://www.numpy.org/>`_ >= 1.6.0
* `SciPy <http://www.scipy.org/>`_ >= 0.10.0
* `Matplotlib <http://matplotlib.org/>`_ >= 1.1
* `Pandas <http://pandas.pydata.org/>`_ >= 0.10.0
* `Astropy <http://www.astropy.org/>`_ >= 0.3.0
* `Suds-jurko <https://bitbucket.org/jurko/suds>`_
* `Beautiful Soup 4 <http://www.crummy.com/software/BeautifulSoup/>`_
* `Glymur <https://github.com/quintusdias/glymur>`_

To import SunPy core (including sun, time, util) you need the following
dependencies:

* `NumPy <http://www.numpy.org/>`_
* `SciPy <http://www.scipy.org/>`_

Any submodules with extra requirements are listed below:

**instr**

* astropy.io.fits
* matplotlib

**io**

* astropy.io.fits, *or*
* glymur

**lightcurve**

* pandas
* astropy.io.fits
* matplotlib

**map**

* astropy.io.fits
* matplotlib

**net**

* suds
* beautifulsoup4

**visualization**

* matplotlib

Installing Scientific Python
----------------------------
For instructions on setting up the scientific Python environment which is
required by SunPy, choose your OS from the list below. When you are done,
return to this page and follow the instructions below to install SunPy and its
requirements.

.. toctree::
   :maxdepth: 1

   mac
   linux
   win

Anaconda
^^^^^^^^
Alternatively, available for all platforms is Anaconda,
a scientific Python distribution that is available free of charge
from `https://store.continuum.io/cshop/anaconda/ <https://store.continuum.io/cshop/anaconda/>`_.
It comes with a complete build environment so you will not need to worry about
installing a compiler or likewise.

Head to the download page, install it and you are set. It can also be installed
into your user directory without administrator permissions; be sure to use
the Anaconda versions of the commands in the following listings if you have
multiple Python environments installed.

Installing Python Modules
-------------------------

Making sure to have followed one of the guides previously. You will have 
installed on your system: Python, Numpy, Scipy, Matplotlib, git, pip and Qt. 
There are few more remaining dependencies which must be installed using pip. 
Depending on which system you are using some of these may already be installed
but it does not hurt to upgrade: ::

 pip install --upgrade distribute
 pip install --upgrade pyfits
 pip install --upgrade suds
 pip install --upgrade pandas
 pip install --upgrade beautifulsoup4

We also recommend you use `IPython <http://ipython.org/>`_ (a great python 
environment). This can also be installed using pip: ::

 pip install --upgrade ipython
 
All done with the SunPy Python prerequisites. You are now ready to install 
SunPy itself.

Installing SunPy
----------------
There are two different versions of SunPy.  Currently, we recommend that users
install the latest stable release.  However, the bleeding edge version on our
GitHub page is quite easy to install and you do get the latest and greatest but
it may have bugs.  Depending on your setup, you may need to preface each of the
``pip ...`` commands with ``sudo pip ...``.

Grabbing the stable SunPy
^^^^^^^^^^^^^^^^^^^^^^^^^
It is as simple as this: ::

    pip install sunpy

If you are upgrading the package: ::

    pip install --upgrade sunpy

For those who like to download the source.
You have a range of download locations.

PyPi: `Download <https://pypi.python.org/packages/source/s/sunpy/sunpy-0.3.1.tar.gz>`_

GitHub (tar.gz): `Download <https://github.com/sunpy/sunpy/tarball/0.3>`__ 

GitHub (zip): `Download <https://github.com/sunpy/sunpy/zipball/0.3>`__ 

Then you can use: ::

    pip install ./<path to download>/SunPyDownload.file_extension

In some cases you may need the ``--no-deps`` flag if pip is trying to upgrade
dependencies such as SciPy and Matplotlib that are difficult to build from
source and the likely errors will abort the upgrade.  That's it folks!

Grabbing the latest SunPy
^^^^^^^^^^^^^^^^^^^^^^^^^
If you do fancy the latest version of SunPy then you will want to `cd` into 
directory you want to have SunPy located and now: ::

 git clone git@github.com:sunpy/sunpy.git

or using HTTP: ::

 git clone http://github.com/sunpy/sunpy.git

This will download the SunPy repository and create a `sunpy` folder at the
current location.  With time, updates will happen and you can update your local
copy by: ::

 git pull upstream master

Finally, to install or upgrade the local version of SunPy you can run: ::

 python setup.py install --upgrade

Testing your installation
-------------------------

Now you can test your installation. Open a new Python shell and type these 
commands: ::

>>> import sunpy
>>> from sunpy.map import Map
>>> Map(sunpy.AIA_171_IMAGE).peek()

If all goes well you should see an AIA 171 image on your screen.

Contributing to SunPy
---------------------
If you are considering contributing to the development of SunPy and please do
consider it.  Please see :ref:`dev-reference-label` as it will explain
everything required to start contributing.

