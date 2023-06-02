.. _sunpy-how-to-create-coordinate-objects:

*************************
Create Coordinate Objects
*************************

The easiest interface to work with coordinates is through the `~astropy.coordinates.SkyCoord` class:

.. code-block:: python

  >>> import astropy.units as u
  >>> from astropy.coordinates import SkyCoord
  >>> from sunpy.coordinates import frames

  >>> coord = SkyCoord(-100*u.arcsec, 500*u.arcsec, frame=frames.Helioprojective)
  >>> coord = SkyCoord(x=-72241.0*u.km, y=361206.1*u.km, z=589951.4*u.km, frame=frames.Heliocentric)
  >>> coord = SkyCoord(70*u.deg, -30*u.deg, frame=frames.HeliographicStonyhurst)
  >>> coord
  <SkyCoord (HeliographicStonyhurst: obstime=None, rsun=695700.0 km): (lon, lat) in deg
      (70., -30.)>

It is also possible to use strings to specify the frame but in that case make sure to explicitly import `sunpy.coordinates` as it registers sunpy's coordinate frames the Astropy coordinates framework:

.. code-block:: python

    >>> import sunpy.coordinates

    >>> coord = SkyCoord(-100*u.arcsec, 500*u.arcsec, frame='helioprojective', observer='earth')
    >>> coord
    <SkyCoord (Helioprojective: obstime=None, rsun=695700.0 km, observer=earth): (Tx, Ty) in arcsec
        (-100., 500.)>

`~astropy.coordinates.SkyCoord` and all coordinate frames support array coordinates.
These work the same as single-value coordinates, but they store multiple coordinates in a single object.
When you're going to apply the same operation to many different coordinates, this is a better choice than a list of `~astropy.coordinates.SkyCoord` objects, because it will be *much* faster than applying the operation to each `~astropy.coordinates.SkyCoord` in a ``for`` loop:

.. code-block:: python

    >>> coord = SkyCoord([-500, 400]*u.arcsec, [100, 200]*u.arcsec, frame=frames.Helioprojective)
    >>> coord
    <SkyCoord (Helioprojective: obstime=None, rsun=695700.0 km, observer=None): (Tx, Ty) in arcsec
        [(-500.,  100.), ( 400.,  200.)]>
    >>> coord[0]
    <SkyCoord (Helioprojective: obstime=None, rsun=695700.0 km, observer=None): (Tx, Ty) in arcsec
        (-500.,  100.)>
