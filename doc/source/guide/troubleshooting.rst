.. _troubleshooting-faq:

***************
Troubleshooting
***************

.. contents::
   :backlinks: none

.. _sunpy-version:

Obtaining sunpy version
============================

To find out your sunpy version number, import it and print the
``__version__`` attribute::

    >>> import sunpy
    >>> matplotlib.__version__
    '0.1'


.. _locating-sunpy-install:

:file:`sunpy` install location
===================================

You can find what directory sunpy is installed in by importing it
and printing the ``__file__`` attribute::

    >>> import sunpy
    >>> sunpy.__file__
    '/home/jdhunter/dev/lib64/python2.5/site-packages/matplotlib/__init__.pyc'

.. _locating-matplotlib-config-dir:

:file:`.sunpy` directory location
======================================

Each user has a :file:`.sunpy/` directory which may contain a
:ref:`sunpyrc <customizing-with-sunpyrc-files>` file. To locate your :file:`.sunpy/`
directory, use :func:`matplotlib.get_configdir`::

    >>> import sunpy as sun
    >>> sun.get_configdir()
    '/home/moon/.matplotlib'

On unix-like systems, this directory is generally located in your
:envvar:`HOME` directory.  On windows, it is in your documents and
settings directory by default::

    >>> import sunpy
    >>> sunpy.get_configdir()
    'C:\\Documents and Settings\\jdhunter\\.sunpy'

If you would like to use a different configuration directory, you can
do so by specifying the location in your :envvar:`SUNPY_CONFIGDIR`
environment variable -- see
:ref:`setting-linux-osx-environment-variables`.


.. _reporting-problems:

Report a problem
================

If you are having a problem with matplotlib, search the mailing
lists first: it is possible that someone else has already run into
your problem.

If not, please provide the following information in your e-mail to the
`mailing list
<http://lists.sourceforge.net/mailman/listinfo/matplotlib-users>`_:

  * your operating system; (Linux/UNIX users: post the output of ``uname -a``)

  * sunpy version::

        python -c `import matplotlib; print matplotlib.__version__`

  * where you obtained sunpy.

  * any customizations to your ``sunpyrc`` file (see
    :ref:`customizing-sunpyrc`).

  * if the problem is reproducible, please try to provide a *minimal*,
    standalone Python script that demonstrates the problem.  This is
    *the* critical step.  If you can't post a piece of code that we
    can run and reproduce your error, the chances of getting help are
    significantly diminished.  Very often, the mere act of trying to
    minimize your code to the smallest bit that produces the error
    will help you find a bug in *your* code that is causing the
    problem.


You will likely get a faster response writing to the mailing list than
filing a bug in the bug tracker.  Most developers check the bug
tracker only periodically.  If your problem has been determined to be
a bug and can not be quickly solved, you may be asked to file a bug in
the tracker so the issue doesn't get lost.