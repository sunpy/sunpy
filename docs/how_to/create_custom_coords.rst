.. _sunpy-how-to-create-custom-coordinate-objects:

*************************
Create Coordinate Objects
*************************

It is possible to specify an observer other than Earth by providing a `SkyCoord` object representing the observer's location. This can be useful when working with data from instruments located at specific points in space, such as spacecraft.

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord, CartesianRepresentation
    >>> from sunpy.coordinates import frames
    >>> import astropy.units as u

Define the observer's position in Heliocentric coordinates

.. code-block:: python

    >>> obserevr_coord = SkyCoord(-100*u.arcsec, 500*u.arcsec, frame=frames.Helioprojective)
    >>> obserevr_coord= SkyCoord(x=-72241.0*u.km, y=361206.1*u.km, z=589951.4*u.km, frame=frames.Heliocentric)
    >>> obserevr_coord = SkyCoord(70*u.deg, -30*u.deg, frame=frames.HeliographicStonyhurst)
    >>> obserevr_coord
    <SkyCoord (HeliographicStonyhurst: obstime=None, rsun=695700.0 km): (lon, lat) in deg
        (70., -30.)>

Create a coordinate in Helioprojective with the custom observer

.. code-block:: python

    >>> coord = SkyCoord(100*u.arcsec, 200*u.arcsec, frame=frames.Helioprojective, observer=observer_coord)
    >>> coord
    <SkyCoord (Helioprojective: obstime=None, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=None, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, km)
        (70., -30., 695700.)>): (Tx, Ty) in arcsec
        (100., 200.)>

In this example, we define the observer's position at -100*u.arcsec along the x-axis and 500*u.arcsec along y-axis in Heliocentric coordinates. This observer coordinate is then used to create a `SkyCoord` object in the `Helioprojective` frame.
