.. _customizing-sunpy:

*****************
Customizing sunpy
*****************

.. _customizing-with-sunpyrc-files:

The :file:`sunpyrc` file
========================

The sunpy core package uses a :file:`sunpyrc` configuration file to customize
certain properties. You can control a number of key features of sunpy such as
where your data will download to. sunpy looks for the ``sunpyrc`` file
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

Do not edit the default file (the first in the "FILES USED:" list above) directly as every time you install or update sunpy, this file will be overwritten.

To maintain your personal customizations place a copy of the default "sunpyrc" file into the user configuration path.
To find this path, use:

.. code-block:: python

    >>> from sunpy.extern.appdirs import AppDirs
    >>> AppDirs('sunpy', 'sunpy').user_config_dir  # doctest: +SKIP

You can use `sunpy.util.config.copy_default_config` to write the default config into the correct place.

The user configuration path can also be set using an environment variable ``SUNPY_CONFIGDIR``.

Depending on your system, it may be useful to have a site-wide configuration file.
If it is used, it will be on the "FILES USED:" list below the default file.
To find your system's site configuration path for ``sunpy``, use:

.. code-block:: python

    >>> from sunpy.extern.appdirs import AppDirs
    >>> AppDirs('sunpy', 'sunpy').site_config_dir  # doctest: +SKIP

In Unix, the site and user configuration paths follow the `XDG specifications <https://standards.freedesktop.org/basedir-spec/basedir-spec-latest.html>`__.

The site configuration is applied before your personal user configuration, thus your configuration file will override the site configuration settings.
For this reason, when you create or edit your personal configuration file, you may want to replicate the site configuration items into your own configuration file, or comment out the items in your configuration file that are set in the site configuration file.

See below for the example config file.

.. _customizing-with-dynamic-settings:

Dynamic settings
===================

You can also dynamically change the default settings in a python script or
interactively from the python shell. All of the settings are stored in a
Python ConfigParser instance called ``sunpy.config``, which is global to
the sunpy package. Settings can be modified directly, for example::

    import sunpy
    sunpy.config.set('downloads', 'download_dir', '/home/user/Downloads')


.. sunpyrc-sample:

A sample sunpyrc file
--------------------------------------------------------------------

.. only:: html

    `(download) <../_static/sunpyrc>`__

.. literalinclude:: ../../sunpy/data/sunpyrc
