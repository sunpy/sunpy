.. _troubleshooting-faq:

***************
Troubleshooting
***************

.. contents::
   :backlinks: none

.. _sunpy-version:

Obtaining sunpy version
=======================

To find out your sunpy version number, import it and print the
``__version__`` attribute::

    >>> import sunpy   # doctest: +SKIP
    >>> sunpy.__version__   # doctest: +SKIP

.. _locating-sunpy-install:

System Info
===========

To quickly collect information on your system, you can use our convenience function
``system_info`` which you can run through: ::

    >>> import sunpy   # doctest: +SKIP
    >>> sunpy.util.system_info()   # doctest: +SKIP

The output should look something like: ::

    ==========================================================
     SunPy Installation Information

     Sunday, 18. November 2012 11:06PM UT
    ==========================================================

    ###########
     General
    ###########
    OS: Mac OS X 10.8.2 (i386)
    Python: 2.7.3 (64bit)

    ####################
     Required libraries
    ####################
    SunPy: 0.1
    NumPy: 1.6.2
    SciPy: 0.10.1
    Matplotlib: 1.2.x
    PyFITS: 3.0.8
    pandas: 0.8.1

    #######################
     Recommended libraries
    #######################
    beautifulsoup4: 4.1.1
    PyQt: 4.9.4
    SUDS: 0.4'

This information is especially useful if you are running into a bug and need help.

:file:`sunpy` install location
===================================

You can find what directory sunpy is installed in by importing it
and printing the ``__file__`` attribute::

    >>> import sunpy   # doctest: +SKIP
    >>> sunpy.__file__   # doctest: +SKIP

.. _locating-matplotlib-config-dir:

:file:`.sunpy` directory location
=================================

Each user should have a :file:`.sunpy/` directory which should contain a
:ref:`sunpyrc <customizing-with-sunpyrc-files>` file. To locate your :file:`.sunpy/`
directory, use :func:`sunpy.print_config`::

    >>> import sunpy as sun   # doctest: +SKIP
    >>> sun.print_config()   # doctest: +SKIP

On unix-like systems, this directory is generally located in your
:envvar:`HOME` directory.  On windows, it is in your documents and
settings directory by default.

If you would like to use a different configuration directory, you can
do so by specifying the location in your :envvar:`SUNPY_CONFIGDIR`
environment variable.
