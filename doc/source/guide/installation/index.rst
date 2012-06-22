============
Installation
============

Below are instructions for installation SunPy and its prerequisites on 
different platforms.

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

Installing Prerequisites
------------------------
For instructions on installing these prereqs on your system, choose your OS from
the list below.

.. toctree::
   :maxdepth: 1
   
   mac
   linux
   win
   
Installing SunPy
----------------
Once you have successfully installed Python and the prereqs as described in the
sections above you are ready to install SunPy itself. There are several 
different methods for installing SunPy. Currently, as SunPy is in active development
we recommend that you download the latest unstable version of SunPy.

.. toctree::
   :maxdepth: 1
   
   stable
   git
   
Testing your installation
-------------------------

Now you can test your installation. Open a new Python shell and type these 
commands: ::

>>> import sunpy
>>> sunpy.make_map(sunpy.AIA_171_IMAGE).show()

If all goes well you should see an AIA 171 image on your screen.
