.. _customizing-sunpy:

*****************
Customizing SunPy
*****************

.. _customizing-with-sunpyrc-files:

The :file:`sunpyrc` file
========================

Sunpy uses a :file:`sunpyrc` configuration file to customize certain
properties. You can control a number of key features of SunPy such as
where your data will download to. SunPy looks for the ``sunpyrc`` file
in a platform specific directory, which you can see the path for by running::

  >>> import sunpy # doctest: +SKIP
  >>> sunpy.print_config()   # doctest: +SKIP

To maintain your own customizations place a copy of the default sunpyrc file
into the *first* path printed above. Do not edit the default file directly as every
time you install or update SunPy, this file will be overwritten.

See below for the example config file.

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

.. only:: html

    `(download) <../_static/sunpyrc>`__

.. literalinclude:: ../../sunpy/data/sunpyrc
