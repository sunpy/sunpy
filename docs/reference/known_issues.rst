.. _known_issues:

************
Known Issues
************

This page documents commonly known issues, issues here is defined broadly and refers to oddities or specifics of how sunpy or the Python ecosystem works that could anyone catch out.

Disagreement between `astropy.time.Time` and SSW for pre-1972 dates
===================================================================

Conversions between time scales for pre-1972 dates may not agree between `astropy.time.Time` and SSW.
Prior to 1972, the length of a Universal Time (UT) second was different from the length of an SI second, and `astropy.time.Time` correctly handles this difference.
However, SSW does not (at least in the commonly used routines), and thus SSW conversions for pre-1972 dates between UT and International Atomic Time (TAI, based on SI seconds) can be incorrect by multiple seconds.
The SSW routine ``utc2tai`` notes that its output is invalid for pre-1972 dates, but other SSW routines are missing such a disclaimer.

`Reference issue <https://github.com/sunpy/sunpy/issues/5500>`__.

`~sunpy.net.jsoc.JSOCClient` time filtering
===========================================

The JSOC API filters on ``T_OBS`` for some series, instead of ``T_REC`` which is what is used most of the time.
There is not anything we can do to fix this on our side.
`The specific section of the SDO users guide explaining this <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/sdoguidese4.html#x9-240004.2.4>`__.

`Reference issue <https://github.com/sunpy/sunpy/issues/5447>`__.
