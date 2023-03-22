.. _troubleshooting-faq:

************************
Troubleshooting and Bugs
************************

.. _sunpy-version:

Obtaining sunpy version
=======================

To find out your sunpy version number, import it and print the ``__version__`` attribute::

    >>> import sunpy   # doctest: +SKIP
    >>> sunpy.__version__   # doctest: +SKIP

.. _locating-sunpy-install:

System Info
===========

To quickly collect information on your system, you can use our convenience function ``system_info`` which you can run through: ::

    >>> import sunpy   # doctest: +SKIP
    >>> sunpy.util.system_info()   # doctest: +SKIP

The output should look something like: ::

    ==========================================================
     sunpy Installation Information

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
    sunpy: 0.1
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

Making use of the logger
========================

For information on configuring and using ``sunpy``\'s logging system, a useful tool for troubleshooting, see :ref:`logger`.

:file:`sunpy` install location
===================================

You can find what directory sunpy is installed in by importing it and printing the ``__file__`` attribute::

    >>> import sunpy   # doctest: +SKIP
    >>> sunpy.__file__   # doctest: +SKIP

.. _locating-matplotlib-config-dir:

:file:`.sunpy` directory location
=================================

Each user should have a :file:`.sunpy/` directory which should contain a :ref:`sunpyrc <customizing-with-sunpyrc-files>` file.
To locate your :file:`.sunpy/` directory, use :func:`sunpy.print_config`::

    >>> import sunpy as sun   # doctest: +SKIP
    >>> sun.print_config()   # doctest: +SKIP

We use `appdirs <https://github.com/ActiveState/appdirs>`__ to work out the location depending on your operating system.

If you would like to use a different configuration directory, you can do so by specifying the location in your ``SUNPY_CONFIGDIR`` environment variable.

.. _reporting-problems:

Reporting Bugs
==============

If you are having a problem with sunpy, search the `mailing list`_ or the github `issue tracker`_.
It is possible that someone else has already run into your problem.

If not, please provide the following information in your e-mail to the `mailing list`_ or to the github `issue tracker`_:

  * your operating system; (Linux/UNIX users: post the output of ``uname -a``)

  * sunpy version::

        >>> import sunpy   # doctest: +SKIP
        >>> sunpy.util.system_info()   # doctest: +SKIP

  * how you obtained sunpy.

  * any customizations to your ``sunpyrc`` file (see :ref:`customizing-sunpy`).

  * Please try to provide a **minimal**, standalone Python script that demonstrates the problem.
    This is **the** critical step.
    If you can't post a piece of code that we can run and reproduce your error, the chances of getting help are significantly diminished.
    Very often, the mere act of trying to minimize your code to the smallest bit that produces the error will help you find a bug in **your** code that is causing the problem.

.. _`mailing list`: https://groups.google.com/forum/#!forum/sunpy
.. _`issue tracker`:  https://github.com/sunpy/sunpy/issues
