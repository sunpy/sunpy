.. _sunpy-coordinates-velocities:

Coordinates with velocity information
*************************************

Velocity information can be added to any coordinate [#differentials]_.
When the coordinate is transformed to a different coordinate frame, the velocity vector will be transformed as appropriate.
Be aware that the transformation framework does not take any velocity information into account when transforming the position vector.

Creating a SkyCoord with velocity
=================================
Velocity information can be added as keyword arguments to `~astropy.coordinates.SkyCoord`.
For SunPy's frames, the names of the velocities components are the names of the position components prepended by "d\_", e.g.,::

    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u
    >>> import sunpy.coordinates

    >>> sc = SkyCoord(lon=10*u.deg, lat=20*u.deg, distance=1*u.AU,
    ...               d_lon=3*u.arcmin/u.week, d_lat=4*u.arcmin/u.d, d_distance=5*u.km/u.min,
    ...               frame='heliocentricinertial', obstime='2021-01-01')
    >>> sc
    <SkyCoord (HeliocentricInertial: obstime=2021-01-01T00:00:00.000): (lon, lat, distance) in (deg, deg, AU)
        (10., 20., 1.)
     (d_lon, d_lat, d_distance) in (arcsec / s, arcsec / s, km / s)
        (0.00029762, 0.00277778, 0.08333333)>

See :ref:`astropy-coordinates-velocities` for ways to add velocity information to existing coordinates.

Querying velocity information
=============================
SunPy has functions to query the positions of planets or other objects (e.g., :func:`~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst`).
For any of these functions, if ``include_velocity=True`` is specified, the returned coordinate will include velocity information, e.g.,::

    >>> from sunpy.coordinates import get_body_heliographic_stonyhurst
    >>> get_body_heliographic_stonyhurst('mercury', '2021-01-01', include_velocity=True)
    <HeliographicStonyhurst Coordinate (obstime=2021-01-01T00:00:00.000): (lon, lat, radius) in (deg, deg, AU)
        (-156.46460438, -1.38836399, 0.43234904)
     (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
        (0.09138133, 0.00720229, -7.2513617)>

Transforming velocities
=======================
The transformation of the velocity vector between two coordinate frames takes into account two effects:

* The change in the direction of the velocity vector due to a change in the orientation of the axes between the two frames.
* The "induced" velocity due to the time dependence of the frames themselves.

Orientation change
------------------
To illustrate the orientation change, let's start with the `~astropy.coordinates.SkyCoord` created at the beginning, which was defined in the `~sunpy.coordinates.frames.HeliocentricInertial` frame.
We transform to Astropy's `~astropy.coordinates.HCRS` frame, which is a different inertial [#inertial]_ frame that is also centered at the Sun::

    >>> sc_hcrs = sc.transform_to('hcrs')
    >>> sc_hcrs
    <SkyCoord (HCRS: obstime=2021-01-01T00:00:00.000): (ra, dec, distance) in (deg, deg, AU)
        (80.95428245, 44.31500877, 1.)
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        (-8828397.4990187, 87659732.13232313, 0.08333338)>

Even though the velocity vectors are oriented very differently in their respective spherical coordinates, their amplitudes are essentially the same::

    >>> sc.velocity.norm()
    <Quantity 2.02654081 km / s>
    >>> sc_hcrs.velocity.norm()
    <Quantity 2.02654083 km / s>

Induced velocity
----------------
To illustrate "induced" velocity, consider the `~sunpy.coordinates.frames.HeliographicStonyhurst` frame, which is defined such that the Earth is always at zero degrees longitude.
That is, this frame rotates around the Sun over time to "follow" the Earth.
Accordingly, the longitude component of Earth's velocity vector will be negligible in this frame::

    >>> from sunpy.coordinates import get_earth
    >>> earth = get_earth('2021-01-01', include_velocity=True)
    >>> earth
    <SkyCoord (HeliographicStonyhurst: obstime=2021-01-01T00:00:00.000): (lon, lat, radius) in (deg, deg, AU)
        (0., -3.02983361, 0.98326486)
     (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
        (1.82278759e-09, -0.00487486, -0.01720926)>

Transforming this coordinate to the `~sunpy.coordinates.frames.HeliocentricInertial` frame, which does not rotate over time, confirms that the Earth is moving in inertial space at the expected ~1 degree/day in heliographic longitude::

    >>> earth.heliocentricinertial
    <SkyCoord (HeliocentricInertial: obstime=2021-01-01T00:00:00.000): (lon, lat, distance) in (deg, deg, AU)
        (24.55623543, -3.02983361, 0.98326486)
     (d_lon, d_lat, d_distance) in (arcsec / s, arcsec / s, km / s)
        (0.0422321, -0.00487486, -0.01720925)>
    >>> earth.heliocentricinertial.d_lon.to('deg/d')
    <Quantity 1.01357048 deg / d>

Transforming over time
======================
As the transformation framework is currently implemented, transforming between frames with different values of ``obstime`` takes into account any time dependency for the definitions of the frames, but does *not* incorporate any notion of the coordinate itself moving in inertial space.
This behavior does not change even if there is velocity information attached to the coordinate.
For example, if we take the same coordinate created earlier for Earth, and transform it to one day later::

    >>> from sunpy.coordinates import HeliographicStonyhurst
    >>> earth.transform_to(HeliographicStonyhurst(obstime=earth.obstime + 1*u.day))
    <SkyCoord (HeliographicStonyhurst: obstime=2021-01-02T00:00:00.000): (lon, lat, radius) in (deg, deg, AU)
        (-1.01325847, -3.02987313, 0.98326043)
     (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
        (-1.15705401e-05, -0.00487486, -0.01698849)>

Note that the location of the Earth in the new frame is ~-1 degree in longitude, as opposed to zero degrees.
That is, this coordinate represents the location of Earth on 2021 January 1 using axes that are defined using the location of Earth on 2021 January 2.

Footnotes
=========

.. [#differentials] Differentials of position with respect to units other than time are also possible, but are not currently well supported.
.. [#inertial] While `~sunpy.coordinates.frames.HeliocentricInertial` and `~astropy.coordinates.HCRS` have fixed axes directions, strictly speaking the small motion of the origin (the Sun) will induce a translational velocity relative to `~astropy.coordinates.ICRS`, but that aspect will cancel out in the transformation.
