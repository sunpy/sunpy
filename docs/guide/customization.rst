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

  >>> import sunpy
  >>> sunpy.print_config()  # doctest: +SKIP
  FILES USED:
    ...
  <BLANKLINE>
  CONFIGURATION:
    [general]
    time_format = %Y-%m-%d %H:%M:%S
    working_dir = ...
  <BLANKLINE>
    [downloads]
    download_dir = ...
    remote_data_manager_dir = ...
    cache_expiry = 10
    sample_dir = ...
  <BLANKLINE>
    [database]
    url = sqlite:////...
  <BLANKLINE>
    [logger]
    log_level = INFO
    use_color = True
    log_warnings = True
    log_exceptions = False
    log_to_file = False
    log_file_level = INFO
    log_file_format = %(asctime)s, %(origin)s, %(levelname)s, %(message)s
  <BLANKLINE>

To maintain your own customizations place a copy of the default sunpyrc file
into the *first* path printed above.
You can use `sunpy.util.config.copy_default_config` to write the default config into the correct place.
Do not edit the default file directly as every time you install or update SunPy, this file will be overwritten.

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
