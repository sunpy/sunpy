.. _customizing-sunpy:

**********************
Customizing sunpy
**********************

.. _customizing-with-sunpyrc-files:

The :file:`sunpyrc` file
=============================

sunpy uses :file:`sunprc` configuration files to customize all kinds
of properties. You can control a number of key features of sunpy such as 
where your data will download. sunpy looks for :file:`sunpy` in three
locations, in the following order:

1. :file:`.sunpy/sunprc`, for the user's default customizations.
2. :file:`{INSTALL}/sunpy/data/sunprc`, where :file:`{INSTALL}`
   is something like :file:`/usr/lib/python2.5/site-packages` on Linux, and
   maybe :file:`C:\\Python25\\Lib\\site-packages` on Windows. Every time you
   install sunpy, this file will be overwritten, so if you want your
   customizations to be saved, please move this file to your :file:`.sunpy`
   directory.

To display where the currently active :file:`sunpyrc` file was
loaded from, one can do the following::

  >>> import sunpy
  >>> sunpy.sunpy_configfile()
  '/home/foo/.sunpy/sunpyrc'

See below for an example config file.

.. _customizing-with-dynamic-settings:

Dynamic settings
===================

You can also dynamically change the default settings in a python script or
interactively from the python shell. All of the settings are stored in a
dictionary-like variable called :data:`sunpy.configs`, which is global to
the sunpy package. Configs can be modified directly, for example::

    import sunpy
    sunpy.config['data.save'] = True
    sunpy.config['data.location'] = '/home/foo/desktop/'


.. sunpyrc-sample:

A sample sunpyrc file
--------------------------------------------------------------------

.. htmlonly::

    `(download) <../_static/sunpyrc>`__

.. literalinclude:: ../../../sunpy/data/sunpyrc