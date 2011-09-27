============
Installation
============

Below are instructions for installation SunPy and its prerequisites on 
different platforms.

Installing Prerequisites
------------------------
**Required**

For its basic functioning, SunPy requires several libraries:

* `NumPy <http://numpy.scipy.org/>`__
* `Matplotlib <http://matplotlib.sourceforge.net/>`__
* `PyFITS <http://www.stsci.edu/resources/software_hardware/pyfits>`_

**Optional**

In addition to the required libraries listed above, there are a couple other
optional dependencies which are only needed for certain features.

* `PyQt4 <http://www.riverbankcomputing.co.uk/software/pyqt/download>`__ (SunPy Plotman)
* `Suds <https://fedorahosted.org/suds/>`__ (VSO/HEK support)

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
different methods for installing SunPy. In general, unless you require a more
recent version or plan to contribute to the development of SunPy, it is 
recommended that you use one of the methods for installing the latest stable
version of SunPy.

.. toctree::
   :maxdepth: 1
   
   stable
   git
   
Testing your installation
-------------------------

Now you can test your installation. Open a new Python shell by typing 
:command:`python` in ``Command Prompt``, and type these commands: ::

>>> import sunpy
>>> sunpy.Map(sunpy.AIA_171_IMAGE).show()