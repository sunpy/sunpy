.. _logger:

**************
Logging system
**************

Overview
========

The SunPy logging system makes use of the `~astropy.logger.AstropyLogger`.
Its purpose is to provide users the ability to decide which log and warning messages to show,
to capture them, and to send them to a file.

All messages provided by SunPy use this logging facility which is based
on the Python `logging` module rather than print statements.

Messages can have one of several levels, in increasing order of importance:

* DEBUG: Detailed information, typically of interest only when diagnosing
  problems.

* INFO: A message conveying information about the current task, and
  confirming that things are working as expected

* WARNING: An indication that something unexpected happened, and that user
  action may be required.

* ERROR: indicates a more serious issue where something failed but the task is continuing

* CRITICAL: A serious error, indicating that the program itself may be unable to continue running.

By default, all messages except for DEBUG messages are displayed.

Configuring the logging system
==============================
The default configuration for the logger is determined by the default SunPy
configuration file. To make permanent changes to the logger configuration
see the ``[logger]`` section of the Sunpy configuration
file (:doc:`sunpyrc </guide/customization>`).

If you'd like to control the logger configuration for your current session
first import the logger::

    >>> from sunpy import log

or also by::

    >>> import logging
    >>> log = logging.getLogger('sunpy')

The threshold level for messages can be set with::

    >>> log.setLevel('DEBUG')  #doctest: +SKIP

This will display DEBUG and all messages with that level and above. If you'd like to see the fewest
relevant messages you'd set the logging level to WARNING or above.

For other options such as whether to log to a file or what level of messages the log file should
contain, see the the Sunpy configuration file (:doc:`sunpyrc </guide/customization>`).

Context managers
================
If you'd like to
capture messages as they are generated you can do that with a context manager::

    >>> from sunpy import log
    >>> with log.log_to_list() as log_list:  #doctest: +SKIP
    ...    # your code here  # doctest: +SKIP

Once your code is executed, ``log_list`` will be a Python list containing all of the Sunpy
messages during execution. This does not divert the messages from going to a file or to the screen.
It is also possible to send the messages to a custom file with::

    >>> from sunpy import log
    >>> with log.log_to_file('myfile.log'):  #doctest: +SKIP
    ...     # your code here  #doctest: +SKIP

which will save the messages to a local file called ``myfile.log``.
