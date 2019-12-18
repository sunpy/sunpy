.. _sunpy-coordinates-rotatedsunframe:

Differential rotation using coordinate frames
=============================================

Normally, coordinates refer to a point in inertial space (relative to the barycenter of the solar system).
Transforming to a different observation time does not move the point at all, but rather only updates the coordinate representation as needed for the origin and axis orientations at the new observation time.
However, the Sun both moves translationally in inertial space and rotates about its axis over time.
Thus, a coordinate that points to a solar feature (e.g., the center of the Sun or an active region on its surface) will not continue to point to the same solar feature at other observation times.

To evolve a coordinate in time such that it accounts for the rotational motion of the Sun, one can use the `~sunpy.coordinates.metaframes.RotatedSunFrame` "metaframe" class as described below.
This machinery will take into account the latitude-dependent rotation rate of the solar surface, also known as differential rotation.

In addition, one may want to account for the translational motion of the Sun as well, and that can be achieved by also using the context manager `~sunpy.coordinates.transform_with_sun_center` for desired coordinate transformations.

Basics of the RotatedSunFrame class
-----------------------------------
In a nutshell, the `~sunpy.coordinates.metaframes.RotatedSunFrame` class allows one to specify coordinates in a coordinate frame prior to an amount of solar (differential) rotation being applied.
That is, the coordinate will point to a location in inertial space at some time, but will use a coordinate system at a *different* time to refer to that point, while accounting for the differential rotation between those two times.

Technical note: ``RotatedSunFrame`` is not itself a coordinate frame, but is instead a "metaframe".
A new frame class is created on the fly corresponding to each base coordinate frame class.
This tutorial will refer to these new classes as ``RotatedSun*`` frames.

Creating coordinates
^^^^^^^^^^^^^^^^^^^^

``RotatedSunFrame`` requires two inputs: the base coordinate frame and the duration of solar rotation.
The base coordinate frame needs to be fully specified, which means a defined ``obstime`` and, if relevant, a defined ``observer``.
Note that the ``RotatedSun*`` frame that is created in this example is appropriately named ``RotatedSunHeliographicStonyhurst``::

  >>> from astropy.coordinates import SkyCoord
  >>> import astropy.units as u
  >>> from sunpy.coordinates import RotatedSunFrame
  >>> import sunpy.coordinates.frames as f

  >>> rs_hgs = RotatedSunFrame(base=f.HeliographicStonyhurst(obstime="2001-01-01"), duration=1*u.day)
  >>> rs_hgs
  <RotatedSunHeliographicStonyhurst Frame (base=<HeliographicStonyhurst Frame (obstime=2001-01-01T00:00:00.000)>, duration=1.0 d, rotation_model=howard)>

Once a ``RotatedSun*`` frame is created, it can be used in the same manner as other frames.  Here, we create a `~astropy.coordinates.SkyCoord` using the ``RotatedSun*`` frame::

  >>> SkyCoord(0*u.deg, 0*u.deg, frame=rs_hgs)
  <SkyCoord (RotatedSunHeliographicStonyhurst: base=<HeliographicStonyhurst Frame (obstime=2001-01-01T00:00:00.000)>, duration=1.0 d, rotation_model=howard): (lon, lat, radius) in (deg, deg, km)
      (0., 0., 695700.)>

Instead of explicitly specifying the duration of solar rotation, one can use the keyword argument ``rotated_time``.
The duration will be automatically calculated from the difference between ``rotated_time`` and the ``obstime`` value of the base coordinate frame.
Here, we also include coordinate data in the supplied base coordinate frame::

  >>> rs_hgc = RotatedSunFrame(base=f.HeliographicCarrington(10*u.deg, 20*u.deg, obstime="2020-03-04 00:00"),
  ...                         rotated_time="2020-03-06 12:00")
  >>> rs_hgc
  <RotatedSunHeliographicCarrington Coordinate (base=<HeliographicCarrington Frame (obstime=2020-03-04T00:00:00.000)>, duration=2.5 d, rotation_model=howard): (lon, lat, radius) in (deg, deg, km)
      (10., 20., 695700.)>

A ``RotatedSun*`` frame containing coordinate data can be supplied to ``SkyCoord`` as normal::

  >>> SkyCoord(rs_hgc)
  <SkyCoord (RotatedSunHeliographicCarrington: base=<HeliographicCarrington Frame (obstime=2020-03-04T00:00:00.000)>, duration=2.5 d, rotation_model=howard): (lon, lat, radius) in (deg, deg, km)
      (10., 20., 695700.)>

It is important to understand that the ``RotatedSun*`` coordinate *already* has differential rotation applied.
If one converts the ``RotatedSun*`` coordinate to a "real" coordinate frame, even the base coordinate frame used in the ``RotatedSun*`` frame, one sees that the longitude for the coordinate is different from the initial representation (in this case, ~45 degrees instead of 10 degrees)::

  >>> rs_hgc.transform_to(rs_hgc.base)
  <HeliographicCarrington Coordinate (obstime=2020-03-04T00:00:00.000): (lon, lat, radius) in (deg, deg, km)
      (45.13354448, 20., 695700.)>

The above examples used the default differential-rotation model, but any of the models available through `sunpy.physics.diff_rot` are selectable.
For example, instead of the default ("howard"), one can specify "allen" using the keyword argument ``rotation_model``.
Note the slight difference in the "real" longitude compared to the output above::

  >>> allen = RotatedSunFrame(base=f.HeliographicCarrington(10*u.deg, 20*u.deg, obstime="2020-03-04 00:00"),
  ...                         rotated_time="2020-03-06 12:00", rotation_model="allen")
  >>> allen.transform_to(allen.base)
  <HeliographicCarrington Coordinate (obstime=2020-03-04T00:00:00.000): (lon, lat, radius) in (deg, deg, km)
      (45.22266666, 20., 695700.)>

Transforming coordinate arrays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For another transformation example, we define a meridan with a Carrington longitude of 100 degrees, plus 1 day of differential rotation.
Again, the coordinates are already differentially rotated in inertial space; the ``RotatedSun*`` frame allows one to represent the coordinates in a frame *prior* to the differential rotation::

  >>> meridian = RotatedSunFrame([100]*11*u.deg, range(-75, 90, 15)*u.deg,
  ...                            base=f.HeliographicCarrington(obstime="2001-01-01"),
  ...                            duration=1*u.day)
  >>> meridian
  <RotatedSunHeliographicCarrington Coordinate (base=<HeliographicCarrington Frame (obstime=2001-01-01T00:00:00.000)>, duration=1.0 d, rotation_model=howard): (lon, lat, radius) in (deg, deg, km)
      [(100., -75., 695700.), (100., -60., 695700.), (100., -45., 695700.),
       (100., -30., 695700.), (100., -15., 695700.), (100.,   0., 695700.),
       (100.,  15., 695700.), (100.,  30., 695700.), (100.,  45., 695700.),
       (100.,  60., 695700.), (100.,  75., 695700.)]>

An easy way to "see" the differential rotation is to transform the coordinates to the base coordinate frame.
Note that the points closer to the equator (latitude of 0 degrees) have evolved farther in longitude than the points at high latitudes::

  >>> meridian.transform_to(meridian.base)
  <HeliographicCarrington Coordinate (obstime=2001-01-01T00:00:00.000): (lon, lat, radius) in (deg, deg, km)
      [(110.7550473 , -75., 695700.), (111.70697161, -60., 695700.),
       (112.80904447, -45., 695700.), (113.68216339, -30., 695700.),
       (114.17617983, -15., 695700.), (114.32632838,   0., 695700.),
       (114.17617983,  15., 695700.), (113.68216339,  30., 695700.),
       (112.80904447,  45., 695700.), (111.70697161,  60., 695700.),
       (110.7550473 ,  75., 695700.)]>

In the specific case of `~sunpy.coordinates.frames.HeliographicCarrington`, this frame rotates with the Sun, but in a non-differential manner.
The Carrington longitude approximately follows the rotation of the Sun.
One can transform to the coordinate frame of 1 day in the future to see the difference between Carrington rotation and differential rotation.
Note that equator rotates slightly faster than the Carrington rotation rate (its longitude is now greater than 100 degrees), but most latitudes rotate slower than the Carrington rotation rate::

  >>> meridian.transform_to(f.HeliographicCarrington(obstime="2001-01-02"))
  <HeliographicCarrington Coordinate (obstime=2001-01-02T00:00:00.000): (lon, lat, radius) in (deg, deg, km)
      [( 96.50790846, -7.48914278e+01, 695909.38247254),
       ( 97.49463421, -5.98996150e+01, 696243.38561962),
       ( 98.60552101, -4.49146239e+01, 696540.31843083),
       ( 99.48029415, -2.99354294e+01, 696780.00438298),
       ( 99.97427245, -1.49606149e+01, 696946.17523076),
       (100.12423256,  1.15298509e-02, 697027.56204653),
       ( 99.97381379,  1.49828932e+01, 697018.64842723),
       ( 99.47937935,  2.99554157e+01, 696920.03844654),
       ( 98.60439605,  4.49309606e+01, 696738.41638666),
       ( 97.49378223,  5.99111888e+01, 696486.10075649),
       ( 96.50761602,  7.48974456e+01, 696180.22046675)]>

Be aware that transformations with a change in ``obstime`` will also contend with a translation of the center of the Sun.
Note that the ``radius`` component above is no longer precisely on the surface of the Sun.
For precise transformations of solar features, one should also use the context manager `~sunpy.coordinates.transformations.transform_with_sun_center` to account for the translational motion of the Sun.
Using the context manager, the ``radius`` component stays as the solar radius as desired::

  >>> from sunpy.coordinates import transform_with_sun_center
  >>> with transform_with_sun_center():
  ...     print(meridian.transform_to(f.HeliographicCarrington(obstime="2001-01-02")))
  <HeliographicCarrington Coordinate (obstime=2001-01-02T00:00:00.000): (lon, lat, radius) in (deg, deg, km)
      [( 96.5706473 , -7.50000000e+01, 695700.),
       ( 97.52257161, -6.00000000e+01, 695700.),
       ( 98.62464447, -4.50000000e+01, 695700.),
       ( 99.49776339, -3.00000000e+01, 695700.),
       ( 99.99177983, -1.50000000e+01, 695700.),
       (100.14192838,  4.79380955e-15, 695700.),
       ( 99.99177983,  1.50000000e+01, 695700.),
       ( 99.49776339,  3.00000000e+01, 695700.),
       ( 98.62464447,  4.50000000e+01, 695700.),
       ( 97.52257161,  6.00000000e+01, 695700.),
       ( 96.5706473 ,  7.50000000e+01, 695700.)]>

Transforming multiple durations of rotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another common use case for differential rotation is to track a solar feature over a sequence of time steps.
Let's track an active region that starts at `~sunpy.coordinates.frames.Helioprojective` coordinates (-123 arcsec, 456 arcsec), as seen from Earth, and we will look both backwards and forwards in time.
We specify a range of durations from -5 days to +5 days, stepping at 1-day increments::

  >>> durations = range(-5, 6, 1)*u.day
  >>> ar_start = f.Helioprojective([-123]*len(durations)*u.arcsec, [456]*len(durations)*u.arcsec,
  ...                              obstime=["2001-01-01"]*len(durations), observer="earth")
  >>> ar = RotatedSunFrame(base=ar_start, duration=durations)
  >>> ar
  <RotatedSunHelioprojective Coordinate (base=<Helioprojective Frame (obstime=['2001-01-01 00:00:00.000' '2001-01-01 00:00:00.000'
   '2001-01-01 00:00:00.000' '2001-01-01 00:00:00.000'
   '2001-01-01 00:00:00.000' '2001-01-01 00:00:00.000'
   '2001-01-01 00:00:00.000' '2001-01-01 00:00:00.000'
   '2001-01-01 00:00:00.000' '2001-01-01 00:00:00.000'
   '2001-01-01 00:00:00.000'], rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>)>, duration=[-5. -4. -3. -2. -1.  0.  1.  2.  3.  4.  5.] d, rotation_model=howard): (Tx, Ty) in arcsec
      [(-123., 456.), (-123., 456.), (-123., 456.), (-123., 456.),
       (-123., 456.), (-123., 456.), (-123., 456.), (-123., 456.),
       (-123., 456.), (-123., 456.), (-123., 456.)]>

Let's convert to the base coordinate frame to reveal the motion of the active region over time::

  >>> ar.transform_to(ar.base)
  <Helioprojective Coordinate (obstime=['2001-01-01 00:00:00.000' '2001-01-01 00:00:00.000'
   '2001-01-01 00:00:00.000' '2001-01-01 00:00:00.000'
   '2001-01-01 00:00:00.000' '2001-01-01 00:00:00.000'
   '2001-01-01 00:00:00.000' '2001-01-01 00:00:00.000'
   '2001-01-01 00:00:00.000' '2001-01-01 00:00:00.000'
   '2001-01-01 00:00:00.000'], rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      [(-865.54956344, 418.10284813, 0.98251245),
       (-794.6736101 , 429.25935934, 0.98154904),
       (-676.99949185, 439.15848306, 0.98069504),
       (-519.35479485, 447.21239117, 0.98000079),
       (-330.98303969, 452.94056372, 0.97950733),
       (-123.        , 456.        , 0.97924388),
       (  92.27675962, 456.20707835, 0.97922605),
       ( 302.0813494 , 453.54935963, 0.9794549 ),
       ( 493.98430821, 448.18638939, 0.97991687),
       ( 656.65386199, 440.43943386, 0.98058459),
       ( 780.54121099, 430.77097352, 0.98141858)]>

Be aware that these coordinates are represented in the ``Helioprojective`` coordinates as seen from Earth at the base time.
Since the Earth moves in its orbit around the Sun, one may be more interested in representing these coordinates as they would been seen by an Earth observer at each time step.
Note that the active region moves slightly slower across the disk of the Sun because the Earth orbits in the same direction as the Sun rotates, thus reducing the apparent rotation of the Sun::

  >>> earth_hpc = f.Helioprojective(obstime=ar.base.obstime + durations, observer="earth")
  >>> ar.transform_to(earth_hpc)
  <Helioprojective Coordinate (obstime=['2000-12-27 00:00:00.000' '2000-12-28 00:00:00.000'
   '2000-12-29 00:00:00.000' '2000-12-30 00:00:00.000'
   '2000-12-31 00:00:00.000' '2001-01-01 00:00:00.000'
   '2001-01-02 00:00:00.000' '2001-01-03 00:00:00.000'
   '2001-01-04 00:00:00.000' '2001-01-05 00:00:00.000'
   '2001-01-06 00:00:00.000'], rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
      [(-847.35767803, 419.05107437, 0.9822281 ),
       (-764.81081635, 428.25724753, 0.9813374 ),
       (-644.29157717, 437.09454986, 0.98056026),
       (-491.84018851, 445.01082595, 0.97993388),
       (-315.11434361, 451.4724754 , 0.97948809),
       (-123.        , 456.        , 0.97924388),
       (  74.8471119 , 458.20167789, 0.97921235),
       ( 268.50338021, 457.80291945, 0.97939415),
       ( 448.29323287, 454.66911128, 0.97977955),
       ( 605.28861971, 448.82020568, 0.98034895),
       ( 731.76302454, 440.43591752, 0.98107395)]>

Transforming into RotatedSun frames
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

So far, all of the examples show transformations with the ``RotatedSun*`` frame as the starting frame.
The ``RotatedSun*`` frame can also be the destination frame, which can be more intuitive in some situations and even necessary in some others (due to API limitations).
Let's use a coordinate from earlier, which represents the coordinate in a "real" coordinate frame::

  >>> coord = rs_hgc.transform_to(rs_hgc.base)
  >>> coord
  <HeliographicCarrington Coordinate (obstime=2020-03-04T00:00:00.000): (lon, lat, radius) in (deg, deg, km)
      (45.13354448, 20., 695700.)>

If we create a ``RotatedSun*`` frame for a different base time, we can represent that same point using coordinates prior to differential rotation::

  >>> rs_frame = RotatedSunFrame(base=f.HeliographicCarrington(obstime=coord.obstime), rotated_time="2020-03-06 12:00")
  >>> rs_frame
  <RotatedSunHeliographicCarrington Frame (base=<HeliographicCarrington Frame (obstime=2020-03-04T00:00:00.000)>, duration=2.5 d, rotation_model=howard)>

  >>> new_coord = coord.transform_to(rs_frame)
  >>> new_coord
  <RotatedSunHeliographicCarrington Coordinate (base=<HeliographicCarrington Frame (obstime=2020-03-04T00:00:00.000)>, duration=2.5 d, rotation_model=howard): (lon, lat, radius) in (deg, deg, km)
      (10., 20., 695700.)>

There coordinates are stored in the ``RotatedSun*`` frame, but it can be useful to "pop off" this extra layer and retain only the coordinate representation in the base coordinate frame.
There is a convenience method called `~sunpy.coordinates.metaframes.RotatedSunFrame.as_base()` to do exactly that.
Be aware the resulting coordinate does *not* point to the same location in inertial space, despite the superficial similarity.
Essentially, the component values have been copied from one coordinate frame to a different coordinate frame, and thus this is not merely a transformation between coordinate frames::

  >>> new_coord.as_base()
  <HeliographicCarrington Coordinate (obstime=2020-03-04T00:00:00.000): (lon, lat, radius) in (deg, deg, km)
      (10., 20., 695700.)>

Advanced uses of RotatedSunFrame
--------------------------------

Here are advanced examples of how to use ``RotatedSunFrame``.

Plotting differentially rotated coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Coordinates in a ``RotatedSun*`` frame can be plotted onto a map in the same manner as a normal coordinate.
Here, we plot the differentially rotated locations of a solar feature from -5 days to +5 days.

.. plot::
   :include-source:

   from astropy.coordinates import SkyCoord
   import astropy.units as u
   import matplotlib.pyplot as plt

   from sunpy.coordinates import RotatedSunFrame
   from sunpy.data.sample import AIA_171_IMAGE
   from sunpy.map import Map

   m = Map(AIA_171_IMAGE)

   fig = plt.figure()
   ax = plt.subplot(projection=m)
   m.plot(clip_interval=(1., 99.95)*u.percent)

   durations = range(-5, 6, 1)*u.day
   point = SkyCoord([187]*len(durations)*u.arcsec, [283]*len(durations)*u.arcsec, frame=m.coordinate_frame)
   diffrot_point = RotatedSunFrame(base=point, duration=durations)

   middle = len(durations) // 2
   ax.plot_coord(diffrot_point[middle], 'ro', fillstyle='none')
   ax.plot_coord(diffrot_point[:middle], 'bo', fillstyle='none')
   ax.plot_coord(diffrot_point[middle+1:], 'bo', fillstyle='none')

   plt.show()

Grid overlay on a Map
^^^^^^^^^^^^^^^^^^^^^

We can use the ``RotatedSun*`` frame to straightforwardly plot grid overlays.
Here, we plot normal ``HeliographicStonyhurst`` longitude lines in white and ``RotatedSunHeliographicStonyhurst`` longitude lines (with 27 days of differential rotation) in blue.

.. plot::
   :include-source:

   import astropy.units as u
   import matplotlib.pyplot as plt

   from sunpy.coordinates import RotatedSunFrame, HeliographicStonyhurst
   from sunpy.data.sample import AIA_171_IMAGE
   from sunpy.map import Map

   m = Map(AIA_171_IMAGE)

   fig = plt.figure()
   ax = plt.subplot(projection=m)
   m.plot(clip_interval=(1., 99.95)*u.percent)

   overlay = ax.get_coords_overlay('heliographic_stonyhurst')
   overlay[0].set_ticks(spacing=15. * u.deg)
   overlay[1].set_ticks(spacing=90. * u.deg)
   overlay.grid(ls='-', color='white')

   rs_hgs = RotatedSunFrame(base=HeliographicStonyhurst(obstime=m.date), duration=27*u.day)
   overlay = ax.get_coords_overlay(rs_hgs)
   overlay[0].set_ticks(spacing=15. * u.deg)
   overlay[1].set_ticks(spacing=90. * u.deg)
   overlay.grid(ls='-', color='blue')

   plt.show()

Be aware that the distorted ``RotatedSunHeliographicStonyhurst`` longitude lines are plotted on the original coordinate frame, and thus use the ``Helioprojective`` coordinate frame for the location of AIA at the original observation time (as opposed to 27 days later).

Differentially rotating a Map
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::
   This example requires `reproject` 0.6 or later to be installed.

Here is an example of reprojecting a map to another time with an observer at a different location.
To achieve the reprojection, one needs to transform from the original frame to the ``RotatedSun*`` version of the output frame with ``rotated_time`` set to the time of the original frame.
When reprojecting across a change in time, it is generally advisable to use the `~sunpy.coordinates.transform_with_sun_center` context manager.

.. plot::
   :include-source:

   from astropy.coordinates import SkyCoord
   import astropy.units as u
   from astropy.wcs import WCS
   import matplotlib.pyplot as plt
   from reproject import reproject_interp

   from sunpy.coordinates import Helioprojective, RotatedSunFrame, HeliographicStonyhurst, transform_with_sun_center
   from sunpy.data.sample import AIA_171_IMAGE
   from sunpy.map import Map, make_fitswcs_header

   # Load the map and define the output time with +5 days of differential rotation
   m = Map(AIA_171_IMAGE)
   in_time = m.date
   out_time = in_time + 5*u.day

   # Define the output real frame with an associated new observer (at Earth)
   out_frame = Helioprojective(observer='earth', obstime=out_time)

   # Define the output RotatedSunFrame (its `duration` will be -5 days)
   rot_frame = RotatedSunFrame(base=out_frame, rotated_time=in_time)

   # Construct a WCS object for the output map with the output RotatedSunFrame specified
   out_shape = m.data.shape
   out_center = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=out_frame)
   header = make_fitswcs_header(out_shape, out_center, scale=u.Quantity(m.scale))
   out_wcs = WCS(header)
   out_wcs.coordinate_frame = rot_frame

   # Reproject the map from the input frame to the output RotatedSunFrame while following Sun center
   with transform_with_sun_center():
       arr, _ = reproject_interp(m, out_wcs, out_shape)

   # Create the output map and preserve the original map's plot settings
   out_warp = Map(arr, out_wcs)
   out_warp.plot_settings = m.plot_settings

   # Plot the output map next to the original map
   fig = plt.figure(figsize=(12, 4))

   ax1 = fig.add_subplot(1, 2, 1, projection=m)
   m.plot(vmin=0, vmax=20000, title='Original map')
   plt.colorbar()

   ax2 = fig.add_subplot(1, 2, 2, projection=out_warp)
   out_warp.plot(vmin=0, vmax=20000,
                 title=f"Reprojected to an Earth observer {(out_time - in_time).to('day')} later")
   plt.colorbar()

   plt.show()

Instead of the steps to construct a custom WCS object, one can instead use the WCS from an actual ``Map`` object at the desired output time (e.g., the actual AIA observation at the later time).
