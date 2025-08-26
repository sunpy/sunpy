.. _sunpy-how-to-observer-by-coordinate:

*********************************
Specify an observer by coordinate
*********************************

For coordinate frames that include the observer location as part of the definition (e.g., `~sunpy.coordinates.Helioprojective`), there are two main ways to specify the observer location (the ``observer`` frame attribute).
If it makes sense to designate a planetary body as the observer, one can specify the observer to be that body using a string (e.g., ``observer='earth'``), and the observer will be understood to be situated at the center of that body at the time of the coordinate (the ``obstime`` frame attribute).
However, if the precise location on a planetary body is important, or if the observer is elsewhere in the space (e.g., a spacecraft), one can specify that observer location using a `~astropy.coordinates.SkyCoord` object.

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord, CartesianRepresentation
    >>> from sunpy.coordinates import frames
    >>> import astropy.units as u

The `~astropy.coordinates.SkyCoord` for the observer location can come from a variety of sources.
If you are working with a `~sunpy.map.Map` with an observer location stored in the metadata, you can retrieve the observer location using the `~sunpy.map.GenericMap.observer_coordinate` property.
If the observer is a major spacecraft, you might be able to obtain its location using :func:`~sunpy.coordinates.get_horizons_coord`.
Here, we will manually specify an observer location in Stonyhurst heliographic coordinates (`~sunpy.coordinates.HeliographicStonyhurst`):

.. code-block:: python

    >>> observer_coord = SkyCoord(70*u.deg, -30*u.deg, 1*u.AU, frame=frames.HeliographicStonyhurst)
    >>> observer_coord
    <SkyCoord (HeliographicStonyhurst: obstime=None, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (70., -30., 1.)>

Now we can create a coordinate in helioprojective coordinates (`~sunpy.coordinates.Helioprojective`) using the observer location defined above:

.. code-block:: python

    >>> coord = SkyCoord(100*u.arcsec, 200*u.arcsec, frame=frames.Helioprojective, observer=observer_coord)
    >>> coord
    <SkyCoord (Helioprojective: obstime=None, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=None, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (70., -30., 1.)>): (Tx, Ty) in arcsec
        (100., 200.)>

Alternatively, you can first create the observer-dependent frame and then use it to create coordinates:

.. code-block:: python

    >>> frame = frames.Helioprojective(observer=observer_coord)
    >>> coord = SkyCoord(100*u.arcsec, 200*u.arcsec, frame=frame)
    >>> coord
    <SkyCoord (Helioprojective: obstime=None, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=None, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (70., -30., 1.)>): (Tx, Ty) in arcsec
        (100., 200.)>
