.. _dev_logger:

********************************
Logging, Warnings and Exceptions
********************************

Overview
========

**sunpy** makes use of a logging system to deal with messages (see :ref:`logger`).
This provides the users and developers the ability to decide which messages to show, to capture them, and to optionally also send them to a file.
The logger will log all messages issued directly to it but also warnings issued through `warnings.warn` as well as exceptions.

The logger is configured as soon as **sunpy** is imported.
You can access it by importing it explicitly::

    from sunpy import log

Messages can be issued directly to it with the following levels and in the following way::

    log.debug("Detailed information, typically of interest only when diagnosing problems.")

    log.info("A message conveying information about the current task, and confirming that
              things are working as expected.")

Any other level should be escalated into a warning or an exception.

A warning should be given in code if the issue is avoidable and the user code could be modified to eliminate the warning.
This is for example what happens with data files that are missing WCS observer metadata.
Exceptions should be raised only when there is no way for the function to proceed.

Regardless of the type (log message or warning or exception) messages should be incomplete sentences that fully describe the issue and how to fix it if possible.

Issuing Warnings
================

**sunpy** warnings are provided by the `sunpy.util` module.
The primary warning which should be used is `sunpy.util.exceptions.SunpyUserWarning`.
For deprecation use `sunpy.util.exceptions.SunpyDeprecationWarning` or `sunpy.util.exceptions.SunpyPendingDeprecationWarning`.

These three warning types have corresponding functions to raise them::

    >>> from sunpy.util.exceptions import warn_user
    >>> from sunpy.util.exceptions import warn_deprecated
    >>> from sunpy.util.exceptions import warn_metadata

These warning functions must be used to interact correctly with the logging system.
A warning can be issued in the following way::

    >>> from sunpy.util.exceptions import warn_user
    >>> warn_user("You have been warned about something you did not do correctly.")  # doctest: +IGNORE_WARNINGS

Raising Exceptions
==================

Raising errors causes the program to halt, so they should be used only when the function can not continue.

Exceptions are raised simply with::

    >>> raise ValueError("The following error occurred.")  #doctest: +SKIP

For more information on exceptions see the :ref:`python:bltin-exceptions` documentation.
