.. _sunpy-topic-guide-role-of-rsun:

***************************
The role of "rsun" in sunpy
***************************

Unlike a solid ball, the Sun is not defined by a single surface, and as such, the position from which a photon is emitted should best be defined not just with a 2 dimensional coordinate but with the inclusion of a height.
The Solar radius is a common choice for this height, and is defined by the IAU as the distance from the center of the Sun to a layer in the photosphere where the optical depth is 2/3.
However, different wavelengths of light can be emitted from different heights in the solar atmosphere, and so the radius at which a photon is emitted is not always the same as the solar radius.
It is therefore a useful convention to define "rsun" attribute for a map, which is the radius at which the emission is presumed to originate.
This has limitations: data will invariably consist of a range of wavelengths of light and even a single wavelength, especially for coronal emission, comes from a variety of heights
The value of "rsun" is intended as an estimate for this height.

Data providers often include this information for end users via the fits keyword: "RSUN_REF", where this keyword is not available, sunpy assumes the standard value of the Solar radius, as defined by `sunpy.sun.constants`.
Whether "rsun" is set by the data provider or not, we emphasize that this is little more than an approximation.
The exact value used for "rsun" has consequences for working in sunpy with coordinates and transforms, it provides 3 dimensional information about the location of a map pixel.
For example, when loading a fits file with `sunpy.map.Map` a ``map.wcs`` is constructed with a value for "rsun" collected for use in coordinate transformations.
When manipulating these coordinates, idiosyncrasies and unexpected behavior might be encountered, these can be preempted by an understanding of behaviors surrounding "rsun".

Transforming between Helioprojective frames
===========================================

It is possible to convert from a `~sunpy.coordinates.frames.Helioprojective` frame with one observer location to another `~sunpy.coordinates.frames.Helioprojective` frame with a different observer location.
The transformation requires the coordinate to be 3D, so if the initial coordinate is only 2D, the default assumption maps that 2D coordinate to the surface of the Sun, at a radius defined by the "rsun" frame attribute.
The conversion can be performed as follows:

.. code-block:: python

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord

    >>> import sunpy.coordinates
    >>> from sunpy.coordinates import frames

    >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, observer="earth", obstime="2017-07-26", frame=frames.Helioprojective)
    >>> hpc_out = sunpy.coordinates.Helioprojective(observer="venus", obstime="2017-07-26")
    >>> hpc2 = hpc1.transform_to(hpc_out)

An example with two maps, named ``aia`` and ``stereo``:

.. code-block:: python

  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=aia.coordinate_frame)  # doctest: +SKIP
  >>> hpc2 = hpc1.transform_to(stereo.coordinate_frame)  # doctest: +SKIP

Reprojecting between frames with different "rsun"
=================================================

When the values of "rsun" from two wcs instances are different, issues with reprojecting between those frames can be encountered.
:meth:`~sunpy.map.GenericMap.reproject_to` by default, enforces a round-trip behavior, the idea being that you should only trust the reprojection if the 2D coordinates from each observer both resolve to the same 3D point (within a pixel volume).
When the values for "rsun" are different, that criterion fails towards the limb.
In other cases, this furthers results in banding as the criterion fails, then succeeds and then fails again.

Changing the value of "rsun" can fix this issue.
Take the example :ref:`sphx_glr_generated_gallery_map_transformations_reprojection_different_observers.py`.
If we run the reprojection twice, once before fixing the discrepancy in the "rsun_ref" metadata and once after:

.. plot::
    :include-source:
    :context: close-figs

    from matplotlib import pyplot as plt
    from sunpy.map import Map
    from sunpy.data.sample import AIA_193_JUN2012, STEREO_A_195_JUN2012

    plt.rcParams['figure.figsize'] = (16, 8)

    map_aia = Map(AIA_193_JUN2012)
    map_euvi = Map(STEREO_A_195_JUN2012)

    outmap1 = map_euvi.reproject_to(map_aia.wcs)
    map_euvi.meta['rsun_ref'] = map_aia.meta['rsun_ref']
    outmap2 = map_euvi.reproject_to(map_aia.wcs)

    fig = plt.figure()
    ax1 = fig.add_subplot(121, projection=outmap1); outmap1.plot(axes=ax1); plt.title('Without rsun Fix')
    ax2 = fig.add_subplot(122, projection=outmap2); outmap2.plot(axes=ax2); plt.title('With rsun Fix')

    plt.show()

We can see the difference in the appearance of the reprojected maps near the limb.
