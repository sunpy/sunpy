.. _coordinates-sunpy:

***********
Coordinates
***********

This section of the guide introduces how coordinates are represented in sunpy.
sunpy makes use of the `astropy.coordinates` module for this task.

In much the same way as `~astropy.units` are used for representing physical quantities, sunpy uses `astropy.coordinates` to represent points in physical space.
This applies to both points in 3D space and projected coordinates in images.

The astropy coordinates module is primarily used through the `~astropy.coordinates.SkyCoord` class, which also makes use of the astropy units system::

  >>> from astropy.coordinates import SkyCoord
  >>> import astropy.units as u

To enable the use of the solar physics specific frames defined in sunpy we also need to import them::

  >>> from sunpy.coordinates import frames

A SkyCoord object to represent a point on the Sun can then be created::

  >>> coord = SkyCoord(70*u.deg, -30*u.deg, obstime="2017-08-01",
  ...                  frame=frames.HeliographicStonyhurst)
  >>> coord
  <SkyCoord (HeliographicStonyhurst: obstime=2017-08-01T00:00:00.000, rsun=695700.0 km): (lon, lat) in deg
      (70., -30.)>

This `~astropy.coordinates.SkyCoord` object can then be transformed to any other coordinate frame defined either in Astropy or sunpy (see :ref:`sunpy-coordinate-systems` for a list of sunpy frames), for example to transform from the original Stonyhurst frame to a Helioprojective frame::

  >>> coord.transform_to(frames.Helioprojective(observer="earth"))
  <SkyCoord (Helioprojective: obstime=2017-08-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, km)
      (769.96270814, -498.89715922, 1.51668773e+08)>


It is also possible to convert three dimensional positions to astrophysical frames defined in Astropy, for example `~astropy.coordinates.ICRS`.

  >>> coord.transform_to('icrs')
  <SkyCoord (ICRS): (ra, dec, distance) in (deg, deg, km)
    (49.84856512, 0.05394699, 1417743.94689472)>


Observer Location
-----------------

Both `~sunpy.coordinates.frames.Helioprojective` and `~sunpy.coordinates.frames.Heliocentric` frames are defined based on the position of the observer.
Therefore to transform either of these frames to a different frame the location of the observer must be known.
The observer can be specified for a coordinate object using the ``observer`` argument to `~astropy.coordinates.SkyCoord`.
For sunpy to calculate the location of Earth or another solar-system body, it must know the time associated with the coordinate; this is specified with the ``obstime`` argument.

Using the observer location it is possible to convert a coordinate as seen by one observer to a coordinate seen by another::

  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, observer="earth",
  ...                 obstime="2017-07-26",
  ...                 frame=frames.Helioprojective)

  >>> hpc1.transform_to(frames.Helioprojective(observer="venus",
  ...                                          obstime="2017-07-26"))
  <SkyCoord (Helioprojective: obstime=2017-07-26T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'venus'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      (-1285.47497992, 106.20918654, 0.72405937)>


Using Coordinates with Maps
---------------------------

.. plot::
   :include-source:

   sunpy Map uses coordinates to specify locations on the image, and to plot overlays on plots of maps.
   When a Map is created, a coordinate frame is constructed from the header information.
   This can be accessed using ``.coordinate_frame``:

   >>> import sunpy.map
   >>> from sunpy.data.sample import AIA_171_IMAGE   # doctest: +REMOTE_DATA
   >>> m = sunpy.map.Map(AIA_171_IMAGE)  # doctest: +REMOTE_DATA
   >>> m.coordinate_frame  # doctest: +REMOTE_DATA
   <Helioprojective Frame (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
       (-0.00406308, 0.04787238, 1.51846026e+11)>)>

   This can be used when creating a `~astropy.coordinates.SkyCoord` object to set the coordinate system to that image:

   >>> from astropy.coordinates import SkyCoord
   >>> import astropy.units as u
   >>> coord = SkyCoord(100 * u.arcsec, 10*u.arcsec, frame=m.coordinate_frame)  # doctest: +REMOTE_DATA
   >>> coord  # doctest: +REMOTE_DATA
   <SkyCoord (Helioprojective: obstime=2011-06-07T06:33:02.770, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
       (-0.00406308, 0.04787238, 1.51846026e+11)>): (Tx, Ty) in arcsec
       (100., 10.)>

   The `~astropy.coordinates.SkyCoord` object can be converted to a PixelPair object using `~sunpy.map.GenericMap.world_to_pixel`:

   >>> pixel_obj = m.world_to_pixel(coord) # doctest: +REMOTE_DATA
   >>> pixel_obj # doctest: +REMOTE_DATA
   PixelPair(x=<Quantity 551.7680511 pix>, y=<Quantity 515.18266871 pix>)

   This `~astropy.coordinates.SkyCoord` object could also be used to plot a point on top of the map:

   >>> import matplotlib.pyplot as plt
   >>> ax = plt.subplot(projection=m)  # doctest: +SKIP
   >>> m.plot()  # doctest: +SKIP
   <matplotlib.image.AxesImage object at ...>
   >>> _ = ax.plot_coord(coord, 'o')  # doctest: +SKIP

For more information on coordinates see :ref:`sunpy-coordinates` section of the :ref:`reference`.
