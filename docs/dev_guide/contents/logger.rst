.. _dev_logger:

*********************************
Logging, Warnings, and Exceptions
*********************************

Overview
========

SunPy makes use of a logging system to deal with messages (see :ref:`logger`). This provides the users and
developers the ability to decide which messages to show, to capture them, and to optionally also send
them to a file. The logger will log all messages issued directly to it but also warnings issued
through `warnings.warn` as well as exceptions.

The logger is configured as soon as sunpy is imported. You can access it
by importing it explicitly::

    from sunpy import log

Messages can be issued directly to it with the following levels and in the following way::

    log.debug("Detailed information, typically of interest only when diagnosing problems.")

    log.info("A message conveying information about the current task, and confirming that
              things are working as expected.")

    log.warning("An indication that something unexpected happened, or indicative of
                 some problem in the near future (e.g. ‘disk space low’).

                 The software is still working as expected.")
    log.error("Due to a more serious problem, the software has not been able to
               perform some function but the task is still continuing.")

    log.critical("A serious error, indicating that the program itself may be unable to
                  continue running. A real error may soon by issued and the task will fail.")

The difference between logging a warning/error/critical compared to issuing a Python warning or raising
an exception are subtle but important.

Use Python `warnings.warn` in library code if the issue is avoidable and the user code should be
modified to eliminate the warning.

Use ``log.warning()`` if there is likely nothing the user can do about the situation, but the event
should still be noted. An example of this might be if the input data are all zeros. This may be unavoidable or
even by design but you may want to let the user know.

True exceptions (not ``log.error()``) should be raised only when there is no way for the function to proceed.

Regardless of the type (log message or warning or exception) messages should be one or two complete sentences
that fully describe the issue and end with a period.

Issuing Warnings
================
SunPy warnings are provided by the `sunpy.util` module. The primary warning which
should be used is `sunpy.util.exceptions.SunpyUserWarning`. For deprecation use `sunpy.util.exceptions.SunpyDeprecationWarning` or
`sunpy.util.exceptions.SunpyPendingDeprecationWarning`.

These warning classes must be used to interact correctly with the logging system.
A warning can be issued in the following way::

    >>> import warnings
    >>> from sunpy.util import SunpyUserWarning
    >>> warnings.warn("You have been warned about something you did not do correctly.", SunpyUserWarning)  # doctest: +IGNORE_WARNINGS

See the section above for a discussion about the distinction between ``log.warn()`` and :meth:`warnings.warn`.

Raising Exceptions
==================
Raising errors causes the program to halt. Likely the primary error that a sunpy developer will
want to use is

* ValueError: should be raised if the input is not what was expected and could not be used. Functions should not return anything (like None) in that case or issue a warning.

Exceptions are raised simply with::

    >>> raise ValueError("The following error occurred.")  #doctest: +SKIP

For more information on exceptions see the :ref:`python:bltin-exceptions` documentation.
