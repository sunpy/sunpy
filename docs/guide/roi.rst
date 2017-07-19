==================
Region of Interest
==================

.. warning:: This module is under development.

An region of interest (ROI) is an object that contains some basic information about a particular time range in the form of string descriptors. For example, an ROI might denote an interval of troublesome instrument data, such as an encounter with the South Antarctic Anomaly (SAA).

1. Creating an ROI
------------------

You can create an ROI object with the following: ::

    >>> from sunpy.roi import *
    >>> result = roi(times=['2011-02-15 04:34:09','2011-02-15 04:48:21'],description='UV occult.',source='LYRA LYTAF')

This creates an roi called result for the specific time range. Querying the newly created ROI gives the following: ::

    >>> result   # doctest: +NORMALIZE_WHITESPACE
    SunPy Region-of-interest (ROI) object
    -------------------------------------
    Source: 		LYRA LYTAF
    Start time:		2011-02-15T04:34:09
    End time: 		2011-02-15T04:48:21
    Event description:	UV occult.

Check out the code reference for the time range object for more information.
