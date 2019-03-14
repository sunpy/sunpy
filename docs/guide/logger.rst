.. _logger:

**************
Logging system
**************

Overview
========

The Sunpy logging system makes use of the `Astropy logging system <http://docs.astropy.org/en/stable/logging.html>`__.
Its purpose is to provide users the ability to decide which log and warning messages to show,
to capture them, and to send them to a file.

All messages provided by Sunpy use this logging facility which is based
on the Python `~logging` module rather than print statements.

Messages can have one of several levels

* DEBUG: Detailed information, typically of interest only when diagnosing
  problems.

* INFO: A message conveying information about the current task, and
  confirming that things are working as expected

* WARNING: An indication that something unexpected happened, and that user
  action may be required.

* ERROR: indicates a more serious issue, including exceptions

By default, only WARNING and ERROR messages are displayed, and are sent to a
log file located at ``~/.sunpy/sunpy.log``.

The logger
==========
The logger is configured as soon as sunpy is imported. You can access it
by importing it explicitly::

    from sunpy import log

or also by::

    import logging
    log = logging.getLogger('sunpy')

You can use it to issue your own messages or to change its defaults. If you'd like to
capture messages as they are generated you can do that with a context manager::

    with log.log_to_list() as log_list:
        # your code here

Once your code is executed, ``log_list`` will be a Python list containing all of the Sunpy
messages during execution. This does not divert the messages from going to a file or to the screen.
It is also possible to send the messages to a custom file with::

    with log.log_to_file('myfile.log'):
        # your code here

which will save the messages to a local file called ``myfile.log``.

User Configuration
==================
Options for the logger can be set in the ``[logger]`` section
of the Sunpy configuration file (:doc:`Customization </guide/customization>`)::

    [logger]

    # Threshold for the logging messages. Logging messages that are less severe
    # than this level will be ignored. The levels are DEBUG, INFO, WARNING,
    # ERROR
    log_level = INFO

    # Whether to use color for the level names
    use_color = True

    # Whether to log warnings.warn calls
    log_warnings = False

    # Whether to log exceptions before raising them
    log_exceptions = True

    # Whether to always log messages to a log file
    log_to_file = False

    # The file to log messages to
    log_file_path = sunpy.log

    # Threshold for logging messages to log_file_path
    log_file_level = INFO

    # Format for log file entries
    log_file_format = %(asctime)s, %(origin)s, %(levelname)s, %(message)s


