************
Known Issues
************

This page documents commonly known issues, issues here is defined broadly and refers to oddities or specifics of how sunpy or the Python ecosystem works that could anyone catch out.

Differences between ``anytim2tai`` and `astropy.time.Time`
==============================================================

* For dates on and after ``1972-01-01 00:00:00 UTC``, this format is equivalent to that returned by the ``anytim2tai`` routine in SSW.
* For dates between ``1958-01-01 00:00:00 UTC`` and ``1972-01-01 00:00:00 UTC``, ``anytim2tai`` returns a constant value of 9 s difference between UTC and TAI time while `astropy.time.Time` returns a value of 0 s difference between UTC and TAI on ``1958-01-01 00:00:00``, with this difference increasing approximately linearly to 10 s on 1972-01-01 00:00:00, where the two approaches agree.
* UT and TAI were synchronized by definition on ``1958 Jan 1 00:00:00``.
  ``anytim2tai`` should return 0 at that time.
* Even if ``anytim2tai`` is intentionally choosing to ignore the difference in lengths of UT seconds and TAI seconds prior to 1972, it should be returning 10 for pre-1972 times, not 9.
  Its implementation acts as if there was a leap second added on 1972 Jan 1, but there was not.

`Reference issue <https://github.com/sunpy/sunpy/issues/5500>`__.

`~sunpy.net.jsoc.JSOCClient` time filtering
===========================================

The JSOC API filters on ``T_OBS`` for some series, instead of ``T_REC`` which is what is used most of the time.
There is not anything we can do to fix this on our side.
`The specific section of the SDO users guide explaining this <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/sdoguidese4.html#x9-240004.2.4>`__.

`Reference issue <https://github.com/sunpy/sunpy/issues/5447>`__.
