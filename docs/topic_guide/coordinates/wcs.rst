.. _sunpy-topic-guide-coordinates-wcs:

*******************
Coordinates and WCS
*******************

The `sunpy.coordinates` sub-package provides a mapping between FITS-WCS CTYPE convention and the coordinate frames as defined in `sunpy.coordinates`.
This is used via the `astropy.wcs.utils.wcs_to_celestial_frame` function, with which the sunpy frames are registered upon being imported.
This list is used by for example, wcsaxes to convert from `astropy.wcs.WCS` objects to coordinate frames.

The `sunpy.map.GenericMap` class creates `astropy.wcs.WCS` objects as ``amap.wcs``, however, it adds some extra attributes to the `~astropy.wcs.WCS` object to be able to fully specify the coordinate frame.
It adds ``heliographic_observer`` and ``rsun``.

If you want to obtain a un-realized coordinate frame corresponding to a `~sunpy.map.GenericMap` object you can do the following:

.. code-block:: python

    >>> import sunpy.map
    >>> from sunpy.data.sample import AIA_171_IMAGE  # doctest: +REMOTE_DATA

    >>> amap = sunpy.map.Map(AIA_171_IMAGE)  # doctest: +REMOTE_DATA
    >>> amap.observer_coordinate  # doctest: +REMOTE_DATA
    <SkyCoord (HeliographicStonyhurst: obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
        (-0.00406308, 0.04787238, 1.51846026e+11)>

which is equivalent to:

.. code-block:: python

    >>> from astropy.wcs.utils import wcs_to_celestial_frame # doctest: +REMOTE_DATA

    >>> wcs_to_celestial_frame(amap.wcs)  # doctest: +REMOTE_DATA
    <Helioprojective Frame (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
        (-0.00406308, 0.04787238, 1.51846026e+11)>)>
