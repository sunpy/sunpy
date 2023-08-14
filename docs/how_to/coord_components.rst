.. _sunpy-how-to-access-coordinate-components:

****************************
Access coordinate components
****************************

Individual coordinates can be accessed via attributes on the `~astropy.coordinates.SkyCoord` object, but the names of the components of the coordinates can depend on the the frame and the chosen representation (e.g., Cartesian versus spherical).

`~sunpy.coordinates.Helioprojective`
====================================

For the helioprojective frame, the theta_x and theta_y components are accessed  ``Tx`` and ``Ty``, respectively:

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u

    >>> from sunpy.coordinates import frames

    >>> c = SkyCoord(-500*u.arcsec, 100*u.arcsec, frame=frames.Helioprojective)
    >>> c.Tx
    <Longitude -500. arcsec>
    >>> c.Ty
    <Latitude 100. arcsec>

`~sunpy.coordinates.Heliocentric`
=================================

Heliocentric is typically used with Cartesian components:

.. code-block:: python

    >>> c = SkyCoord(-72241.0*u.km, 361206.1*u.km, 589951.4*u.km, frame=frames.Heliocentric)
    >>> c.x
    <Quantity -72241. km>
    >>> c.y
    <Quantity 361206.1 km>
    >>> c.z
    <Quantity 589951.4 km>

`~sunpy.coordinates.HeliographicStonyhurst` and `~sunpy.coordinates.HeliographicCarrington`
===========================================================================================

Both of the heliographic frames have the components of latitude, longitude and radius:

.. code-block:: python

    >>> c = SkyCoord(70*u.deg, -30*u.deg, 1*u.AU, frame=frames.HeliographicStonyhurst)
    >>> c.lat
    <Latitude -30. deg>
    >>> c.lon
    <Longitude 70. deg>
    >>> c.radius
    <Distance 1. AU>

Heliographic Stonyhurst, when used with Cartesian components, is known as Heliocentric Earth Equatorial (HEEQ).
Here's an example of how to use `~sunpy.coordinates.frames.HeliographicStonyhurst` for HEEQ coordinates:

.. code-block:: python

    >>> c = SkyCoord(-72241.0*u.km, 361206.1*u.km, 589951.4*u.km,
    ...              representation_type='cartesian', frame=frames.HeliographicStonyhurst)
    >>> c.x
    <Quantity -72241. km>
    >>> c.y
    <Quantity 361206.1 km>
    >>> c.z
    <Quantity 589951.4 km>
