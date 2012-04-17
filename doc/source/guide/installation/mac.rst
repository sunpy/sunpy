========
Mac OS X
========

Overview
--------
**Required**

For its basic functioning, SunPy requires several libraries:

* `NumPy <http://numpy.scipy.org/>`__
* `SciPy <http://www.scipy.org/>`__
* `Matplotlib <http://matplotlib.sourceforge.net/>`__ (1.0+)
* `PyFITS <http://www.stsci.edu/resources/software_hardware/pyfits>`_
* `Suds <https://fedorahosted.org/suds/>`__
* `pandas <http://pandas.pydata.org/>`_

**Optional**

In addition to the required libraries listed above, if you plan to work with
JPEG 2000 data, you must also install:

* `OpenJPEG <http://www.openjpeg.org/>`__
* `PIL <http://www.pythonware.com/products/pil/>`__

For improved GUI support, you will need to install PyQt4:

* `PyQt4 <http://www.riverbankcomputing.co.uk/software/pyqt/download>`__

Using Macports
--------------
Macports is a free/open source package management system that simplifies the 
installation of other free/open source software referred to as ports. Macports 
installs itself and its ports in /opt/local. Activating a port places links to 
the port software in the expected location on the drive wherever that may be. 
Due to the fact that Mac Os X already comes with its own version of Python 
(which can change without warning during a software update to the OS) installing
SunPy and its dependencies can be a challenge. Macports provides a good solution
to this problem. This tutorial follows the one provided by Thomas Robitaille 
(astrofrog on Github). If you run into any trouble make sure to check out the 
Problem Solving section at the end. This method has been tested on Mac OSX 10.5,
10.6, and 10.7.

Before you start you'll need to install Xcode which is the freely available 
developer environment provided by Apple. You can find it in the Mac App Store.  
Just a warning, it is a rather large download but provides compilers which are 
necessary for Macports. After installing Xcode you can go ahead, download and 
install MacPorts. When that is done open a terminal and update the package index
using the command: ::

  sudo port selfupdate

If you run into problems with the selfupdate, refer to the Problem Solving 
section below (Issue #1 & Issue #2).

You can then install NumPy and Matplotlib using the command: ::

  sudo port install py27-matplotlib py27-numpy

When that is done let's install SciPy: ::

  sudo port install py27-scipy

If you run into an error with the gcc44 build during the scipy installation,
refer to the Problem Solving section below (Issue #3).

Install the last two remaining dependency for SunPy which is pyFits: ::

  sudo port install py27-pyfits py27-suds
 
Another very useful (and recommended tool) is iPython which provides an advanced
Python shell with: ::

  sudo port install py27-ipython

If you are planning on doing any GUI development you'll want to install PyQT 
(a particular tough install without Macports) with the following: ::

  sudo port install py27-pyqt4

Configuration
-------------

Now that the installation is done you'll want to setup a few things. First, 
create a folder called .matplotlib in your home directory (if it does not 
already exist). ::

  mkdir ~/.matplotlib

Copy the default configuration file provided by matplotlib into this directory 
with the following (one-line) command: ::

  cp /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/matplotlib
  /mpl-data/matplotlibrc ~/.matplotlib/

Now edit .matplotlib/matplotlibrc and change the backend from Agg to MacOSX. ::

  cd ~/.matplotlib
  Open matplotlibrc with your favorite editor.
  Find the line that says "backend     : Agg"
  Change to "backend       : MacOSX".

Finally you'll want to set your new installation of python as the default python
(and ipython) to be called with the following command: ::

  sudo port select --set python python27
  sudo port select --set ipython ipython27
 
Validation
----------
You can test if everything worked with the following example. First open ipython
by just typing in the terminal: (You will probably have to completely restart 
the terminal.)  ::

  ipython

Then just cut and paste the following code: ::

  import numpy as np
  import scipy.special as sp
  import matplotlib.pyplot as plt

  x = np.arange(200).astype('float')/10
  y = sp.j0(x)
  plt.plot(x,y)
  plt.title("Bessel Function n = 0")
  plt.show()

This example tests out numpy, scipy, and matplotlib.

Installing other packages
-------------------------
To install packages that are not in macports make sure to use: ::

  python setup.py install --user

(Run this command as written from the directory that contains the package in 
question.)  
This will install the packages in ~/Library/Python/2.7/lib/python/site-packages 
where they will automatically recognized by Python. This will maintain the 
integrity of the the MacPorts file structure. In general, do not install 
anything into /opt/local without using the ports command.

If you would like to use easy_install then remember to set the directory 
manually so that it installs the library into your local directory. Here is an 
example for installing pIDLy: ::

  easy_install --install-dir='~/Library/Python/2.7/lib/python/site-packages' pidly

It is not necessary to use sudo for this command.

Problem Solving
---------------

1) If you installed MacPorts and are getting an error during the selfupdate 
process involving sqlite, try removing the MacPorts directory entirely and 
reinstall.  The directory to remove is /opt/local.  This must be done from the 
terminal.  You may want to make a backup tar file of the directory before 
deleting it.  (Note that there are lots of files.)

2) If during the selfupdate process you get an error with syncing index(es), 
you may be behind a firewall for your rsync port.  To get around this, do the 
following: ::

  cd /opt/local/etc/macports/

Use your favorite editor to open the sources.conf file.
Make the following changes to the file:  ::

  #rsync://rsync.macports.org/release/ports/ [default]
  http://www.macports.org/files/ports.tar.gz [default]

Now *instead* of using sudo port selfupdate, use the following command: ::

  port -d sync

Now move onto the next step (sudo port install py27-matplotlib py27-numpy)...

3) During the installation of scipy, you may run into trouble with building 
gcc44. The following error message may appear: ::

  --->  Building gcc44
  Error: Target org.macports.build returned: shell command failed (see log for details)
  Error: Failed to install gcc44
  Log for gcc44 is at:   /opt/local/var/macports/logs/_opt_local_var_macports_sources_rsync.macports.org
  _release_tarballs_ports_lang_gcc44/gcc44/main.log
  Error: The following dependencies were not installed: gcc44 swig-python bison gsed swig pcre
  Error: Status 1 encountered during processing.
  To report a bug, see <http://guide.macports.org/#project.tickets>

This issue has been noticed by others (https://trac.macports.org/ticket/25713). 
Thankfully there is a simple solution,  just run the following command to clean
up this failed installation: ::

  sudo port clean gcc44

and then run the last command again: ::

  sudo port install py27-scipy

This should now install without any problems. Now move onto the next step (sudo 
port install py27-pyfits)...

Updating
--------

As new versions of matplotlib or scipy are released every once in a while it is necessary to update. Thankfully
macports is built to make this easy. You can do a full upgrade of all of the software that macports installed with 
the following command: ::

  port upgrade installed

and make sure to get yourself a cup of coffee after hitting return as this will probably run for a while. You can also
upgrade individual packages (and their dependencies) with a similar line of code, namely: ::

  port upgrade packagename

You may have to precede those commands with sudo depending on what level of privileges you have on your system. 