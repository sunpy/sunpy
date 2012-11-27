.. _customizing-sunpy:

**********************
Customizing sunpy
**********************

.. _customizing-with-sunpyrc-files:

The :file:`sunpyrc` file
=============================

Sunpy uses :file:`sunpyrc` configuration files to customize certain
properties. You can control a number of key features of SunPy such as 
where your data will download to. SunPy looks for :file:`sunpyrc` in two
locations, in the following order:

1. :file:`.sunpy/sunpyrc`, for the user's default customizations.
2. :file:`{INSTALL}/sunpy/data/sunpyrc`, where :file:`{INSTALL}`
   is something like :file:`/usr/lib/python2.5/site-packages` on Linux, and
   maybe :file:`C:\\Python25\\Lib\\site-packages` on Windows. Every time you
   install SunPy, this file will be overwritten, so if you want your
   customizations to be saved, please move this file to your :file:`.sunpyrc`
   directory.

To display where the currently active :file:`sunpyrc` file was
loaded from, one can do the following::

  >>> import sunpy
  >>> sunpy.print_config()

The default sunpyrc file lives in sunpy/data/. We recommend that you make a copy of this
file and put it where Sunpy will find it. See below for the example config file.

.. _customizing-with-dynamic-settings:

Dynamic settings
===================

You can also dynamically change the default settings in a python script or
interactively from the python shell. All of the settings are stored in a
Python ConfigParser instance called :data:`sunpy.config`, which is global to
the sunpy package. Settings can be modified directly, for example::

    import sunpy
    sunpy.config.set('downloads', 'download_dir', '/home/user/Downloads')


.. sunpyrc-sample:

A sample sunpyrc file
--------------------------------------------------------------------

.. htmlonly::

    `(download) <../_static/sunpyrc>`__

.. literalinclude:: ../../../sunpy/data/sunpyrc