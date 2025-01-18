.. _sunpy-how-to-create-Coordinate-frames-that-depend-upon-Earths-magnetic-dipole:

**************************
Creating Coordinate Frames
**************************

There are 3 coordinate frames that depend upon Earth's magnetic dipole:
:class:`~sunpy.coordinates.Geomagnetic`,
:class:`~sunpy.coordinates.SolarMagnetic`, and
:class:`~sunpy.coordinates.GeocentricSolarMagnetospheric`.

Geomagnetic Coordinate Frame
****************************

A coordinate or frame in the Geomagnetic (MAG) system.

This section demonstrates how to create and use the :class:`~sunpy.coordinates.Geomagnetic` frame.

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord, CartesianRepresentation
    >>> import astropy.units as u
    >>> from sunpy.coordinates.frames import Geomagnetic

.. code-block:: python

    >>> sc = SkyCoord(50 * u.deg, 75 * u.deg,
    ...               obstime="2011/01/05T00:00:50", distance=5 * u.km,
    ...               magnetic_model='igrf13', frame=Geomagnetic)
    >>> sc
    <SkyCoord (Geomagnetic: obstime=2011-01-05T00:00:50.000, magnetic_model=igrf13): (lon, lat, distance) in (deg, deg, km)
        (50., 75., 5.)>

.. code-block:: python

    >>> sc = SkyCoord(CartesianRepresentation(10 * u.km, 1 * u.km, 2 * u.km),
    ...               obstime="2011/01/05T00:00:50", magnetic_model='igrf13', frame=Geomagnetic)
    >>> sc
    <SkyCoord (Geomagnetic: obstime=2011-01-05T00:00:50.000, magnetic_model=igrf13): (lon, lat, distance) in (deg, deg, km)
        (5.71059314, 11.25523973, 10.24695077)>

SolarMagnetic Coordinate Frame
******************************

A coordinate or frame in the Solar Magnetic (SM) system.

This section demonstrates how to create and use the :class:`~sunpy.coordinates.SolarMagnetic` frame.

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord, CartesianRepresentation
    >>> import astropy.units as u
    >>> from sunpy.coordinates.frames import SolarMagnetic

.. code-block:: python

    >>> sc = SkyCoord(50 * u.deg, 75 * u.deg,
    ...               obstime="2011/01/05T00:00:50", distance=5 * u.km,
    ...               magnetic_model='igrf13', frame=SolarMagnetic)
    >>> sc
    <SkyCoord (SolarMagnetic: obstime=2011-01-05T00:00:50.000, magnetic_model=igrf13): (lon, lat, distance) in (deg, deg, km)
        (50., 75., 5.)>

.. code-block:: python

    >>> sc = SkyCoord(CartesianRepresentation(10 * u.km, 1 * u.km, 2 * u.km),
    ...               obstime="2011/01/05T00:00:50", magnetic_model='igrf13', frame=SolarMagnetic)
    >>> sc
    <SkyCoord (SolarMagnetic: obstime=2011-01-05T00:00:50.000, magnetic_model=igrf13): (lon, lat, distance) in (deg, deg, km)
        (5.71059314, 11.25523973, 10.24695077)>

GeocentricSolarMagnetospheric Coordinate Frame
***********************************************

A coordinate or frame in the GeocentricSolarMagnetospheric (GSM) system.

This section demonstrates how to create and use the :class:`~sunpy.coordinates.GeocentricSolarMagnetospheric` frame.

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord, CartesianRepresentation
    >>> import astropy.units as u
    >>> from sunpy.coordinates.frames import GeocentricSolarMagnetospheric

.. code-block:: python

    >>> sc = SkyCoord(50 * u.deg, 75 * u.deg,
    ...               obstime="2011/01/05T00:00:50", distance=5 * u.km,
    ...               magnetic_model='igrf13', frame=GeocentricSolarMagnetospheric)
    >>> sc
    <SkyCoord (GeocentricSolarMagnetospheric: obstime=2011-01-05T00:00:50.000, magnetic_model=igrf13): (lon, lat, distance) in (deg, deg, km)
        (50., 75., 5.)>

.. code-block:: python

    >>> sc = SkyCoord(CartesianRepresentation(10 * u.km, 1 * u.km, 2 * u.km),
    ...               obstime="2011/01/05T00:00:50", magnetic_model='igrf13', frame=GeocentricSolarMagnetospheric)
    >>> sc
    <SkyCoord (GeocentricSolarMagnetospheric: obstime=2011-01-05T00:00:50.000, magnetic_model=igrf13): (lon, lat, distance) in (deg, deg, km)
        (5.71059314, 11.25523973, 10.24695077)>
