.. _sunpy-how-to-transform-between-coordinate-frames:

***********************************
Transform between coordinate frames
***********************************

Both `~astropy.coordinates.SkyCoord` and `~astropy.coordinates.BaseCoordinateFrame` instances have a `~astropy.coordinates.SkyCoord.transform_to` method.
This can be used to transform the frame to any other frame, either implemented in sunpy or in Astropy (see also :ref:`astropy-coordinates-transforming`).
An example of transforming the center of the solar disk to Carrington coordinates is:

.. code-block:: python

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord

    >>> from sunpy.coordinates import frames

    >>> coord = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=frames.Helioprojective, obstime="2017-07-26",
    ...              observer="earth")
    >>> coord
    <SkyCoord (Helioprojective: obstime=2017-07-26T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty) in arcsec
        (0., 0.)>
    >>> coord.transform_to(frames.HeliographicCarrington)
    <SkyCoord (HeliographicCarrington: obstime=2017-07-26T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (lon, lat, radius) in (deg, deg, AU)
        (283.95956776, 5.31701821, 0.00465047)>

It is also possible to transform to any coordinate system implemented in Astropy.
This can be used to find the position of the solar limb in AltAz equatorial coordinates:

.. code-block:: python

    >>> from astropy.coordinates import EarthLocation, AltAz

    >>> time = '2017-07-11 15:00'
    >>> greenbelt = EarthLocation(lat=39.0044*u.deg, lon=-76.8758*u.deg)
    >>> greenbelt_frame = AltAz(obstime=time, location=greenbelt)
    >>> west_limb = SkyCoord(900*u.arcsec, 0*u.arcsec, frame=frames.Helioprojective,
    ...                      observer=greenbelt.get_itrs(greenbelt_frame.obstime), obstime=time)  # doctest: +REMOTE_DATA
    >>> west_limb.transform_to(greenbelt_frame)  # doctest: +REMOTE_DATA
    <SkyCoord (AltAz: obstime=2017-07-11 15:00:00.000, location=(1126916.53031967, -4833386.58391627, 3992696.62211575) m, pressure=0.0 hPa, temperature=0.0 deg_C, relative_humidity=0.0, obswl=1.0 micron): (az, alt, distance) in (deg, deg, m)
        (111.40782056, 57.1660434, 1.51859559e+11)>
