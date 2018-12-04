Sunpy 0.9.5 (2018-11-27)
========================

Bug Fixes
---------

- Timeseries and lightcurve will now respect updated config values for download directory. (`#2844 <https://github.com/sunpy/sunpy/pull/2844>`__)
- Always use _default_wrap_angle rather than hard coding a wrap angle in the init
  of a frame (`#2853 <https://github.com/sunpy/sunpy/pull/2853>`__)
- Ensure imageanimators only slice arrays with integers (`#2856 <https://github.com/sunpy/sunpy/pull/2856>`__)


Sunpy 0.9.4 (2018-11-14)
========================

Features
--------

- Now able to create a `sunpy.map.Map` using an array and a `astropy.wcs.WCS` object. (`#2793 <https://github.com/sunpy/sunpy/pull/2793>`__)


Bug Fixes
---------

- Fix RHESSI obssum file downloading to include the final day in the time range. (`#2714 <https://github.com/sunpy/sunpy/pull/2714>`__)
- User can convert between HPC and HCC coordinates with different observers. This is implemented by automatically transforming the coordinate into HGS and then changing observer, and then transforming back to HCC. (`#2754 <https://github.com/sunpy/sunpy/pull/2754>`__)
- Changed default file type for Helioviewer to prevent decode errors. (`#2771 <https://github.com/sunpy/sunpy/pull/2771>`__)
- Fixed loading in LASCO C3 data. (`#2775 <https://github.com/sunpy/sunpy/pull/2775>`__)
- Increase figure size to avoid cutting off longer colormap names in `sunpy.cm.show_colormaps`. (`#2824 <https://github.com/sunpy/sunpy/pull/2824>`__)
- The sample data directory will no longer be created until files are downloaded
  to it. (`#2836 <https://github.com/sunpy/sunpy/pull/2836>`__)


Improved Documentation
----------------------

- Minor changes to the developer guide regarding sprint labels. (`#2765 <https://github.com/sunpy/sunpy/pull/2765>`__)
- Copyedited and corrected the solar cycles example. (`#2770 <https://github.com/sunpy/sunpy/pull/2770>`__)
- Changed "online" mark to "remote_data" and made formatting of marks consistent. (`#2799 <https://github.com/sunpy/sunpy/pull/2799>`__)
- Add a missing plot to the end of the units and coordinates guide. (`#2813 <https://github.com/sunpy/sunpy/pull/2813>`__)


Trivial/Internal Changes
------------------------

- Miscellaneous fixes to developer docs about building sunpy's documentation. (`#2825 <https://github.com/sunpy/sunpy/pull/2825>`__)
- Changed sunpy.instr.aia.aiaprep to update BITPIX keyword to reflect the float64 dtype. (`#2831 <https://github.com/sunpy/sunpy/pull/2831>`__)
- Fix SunPy Coordinate tests with Astropy 3.1 (`#2838 <https://github.com/sunpy/sunpy/pull/2838>`__)


Sunpy 0.9.3 (2018-09-12)
========================

Bug Fixes
---------

- Correctly import `~astropy.units.allclose` based on astropy version. This means that `sunpy.coordinates` will import when pytest is not installed. (`#2702 <https://github.com/sunpy/sunpy/pull/2702>`__)
- Raise an error when transforming between HPC and HCC frames if the observer is not the same. (`#2725 <https://github.com/sunpy/sunpy/pull/2725>`__)
- Do not attempt to save a FITS header comment for a keyword which is not in the header. This prevents an error on saving some maps after the metadata had been modified but not the comments. (`#2748 <https://github.com/sunpy/sunpy/pull/2748>`__)
- Add support for `HMIMap` objects as input to `sunpy.instr.aia.aiaprep()`. (`#2749 <https://github.com/sunpy/sunpy/pull/2749>`__)


Improved Documentation
----------------------

- Add contribution guidelines for the sunpy example gallery. (`#2682 <https://github.com/sunpy/sunpy/pull/2682>`__)
- Clean up the docstring for `sunpy.physics.differential_rotation.solar_rotate_coordinate` to make the example clearer. (`#2708 <https://github.com/sunpy/sunpy/pull/2708>`__)


Sunpy 0.9.2 (2018-07-27)
========================

Bug Fixes
---------

- Fix the bug which crashes LASCOMap for when 'date-obs' is reformatted agian from a self applied function. (`#2700 <https://github.com/sunpy/sunpy/pull/2700>`__)
- Correctly import `~astropy.units.allclose` based on astropy version. This means that `sunpy.coordinates` will import when pytest is not installed. (`#2702 <https://github.com/sunpy/sunpy/pull/2702>`__)


Sunpy 0.9.1 (2018-07-26)
========================

Features
--------

- MapCube has been renamed to MapSequence.
  MapCube does not actually work on cubes but a sequence of maps.
  Due to this change, MapCube is now deprecated but will not be removed until SunPy 1.0. (`#2603 <https://github.com/sunpy/sunpy/pull/2603>`__)


Bug Fixes
---------

- parse_time now parses numpy.datetime64 correctly (`#2572 <https://github.com/sunpy/sunpy/pull/2572>`__)
- Fix the bug that prevented VSO queries for HMI data from downloading file
  without speicifying `a.Physobs`. (`#2621 <https://github.com/sunpy/sunpy/pull/2621>`__)
- Fix `sunpy.map.mapcube.MapCube.plot`. The code had not been updated to support the changes to the wcsaxes helper functions. (`#2627 <https://github.com/sunpy/sunpy/pull/2627>`__)
- Replace all use of the deprecated ``sunpy.cm.get_cmap`` with ``plt.get_cmap`` to prevent deprecation warnings being raised. (`#2635 <https://github.com/sunpy/sunpy/pull/2635>`__)
- Fix generation of the coordinate transformation graph with Astropy 3.1.dev (`#2636 <https://github.com/sunpy/sunpy/pull/2636>`__)
- Prevent helioviewer from erroring when downloading file to a directory that
  does not exist. It will now create the directory when required. (`#2642 <https://github.com/sunpy/sunpy/pull/2642>`__)
- Fix transformations into/out of Heliographic Stonyhurst frame when
  the coordinate representation is Cartesian. (`#2646 <https://github.com/sunpy/sunpy/pull/2646>`__)
- Support passing Python file objects to `sunpy.io.fits.write`. (`#2688 <https://github.com/sunpy/sunpy/pull/2688>`__)
- Added DRMS to setup.py so sunpy[all] installs it as a dependancy. (`#2693 <https://github.com/sunpy/sunpy/pull/2693>`__)
- Fix eve 0cs timeseries seperator regex to support Python 3.7 (`#2697 <https://github.com/sunpy/sunpy/pull/2697>`__)


Improved Documentation
----------------------

- Organise the gallery into sections based on example type and tidy up a little. (`#2624 <https://github.com/sunpy/sunpy/pull/2624>`__)
- Added gallery example showing the conversion of Helioprojective Coordinates to Altitude/Azimuth Coordinates to and back. (`#2656 <https://github.com/sunpy/sunpy/pull/2656>`__)


Trivial/Internal Changes
------------------------

- Revert the handling of ``quantity_allclose`` now that `astropy/astropy#7252 <https://github.com/astropy/astropy/pull/7252>`__ is merged. This also bumps the minimum astropy 3 version to 3.0.2. (`#2598 <https://github.com/sunpy/sunpy/pull/2598>`__). We still support Astropy 2.
- Sort the ana C source files before building to enable reproducible builds. (`#2637 <https://github.com/sunpy/sunpy/pull/2637>`__)
- We are now using `towncrier <https://github.com/hawkowl/towncrier>`__ to
  generate our changelogs. (`#2644 <https://github.com/sunpy/sunpy/pull/2644>`__)
- Use of ``textwrap`` to keep source code indented when multiline texts is used (`#2671 <https://github.com/sunpy/sunpy/pull/2671>`__)


0.9.0
=====

New Features
------------

- Added TimeUTime class to support utime. [#2409]
- Example for fine-grained use of ticks and grids [#2435]
- Maintiners Workflow Guide [#2411]
- Decorator to append and/or prepend doc strings [#2386]
- Adding `python setup.py test --figure-only` [#2557]
- Fido.fetch now accepts pathlib.Path objects for path attribute.[#2559]
- The `~sunpy.coordinates.HeliographicStonyhurst` coordinate system can now be specified
  using a cartesian system, which is sometimes known as the
  "Heliocentric Earth equatorial" (HEEQ) coordinate system. [#2437]

API Changes
-----------

- ``sunpy.coordinates.representation`` has been removed. Longitude wrapping is now done in the constructor of the frames. [#2431]
- Propagation of ``obstime`` in the coordinate frame transformation has changed, this means in general when transforming directly between frames (not
  ``SkyCoord``) you will have to specify ``obstime`` in more places. [#2461]
- Transforming between Heliographic Stonyhurst and Carrington now requires that ``obstime`` be defined and the same on both the input and output frames. [#2461]
- Removed the figure return from .peek() [#2487]

Bug Fixes
---------

- Improve TimeSeriesBase docstring [#2399]
- Validate that pytest-doctestplus is installed [#2388]
- Fix use of self.wcs in plot in mapbase [#2398]
- Updated docstring with pointer to access EVE data for other levels [#2402]
- Fix broken links and redirections in documentation [#2403]
- Fixes Documentation changes due to NumPy 1.14 [#2404]
- Added docstrings to functions in dowload.py [#2415]
- Clean up database doc [#2414]
- rhessi.py now uses sunpy.io instead of astropy.io [#2416]
- Remove Gamma usage in Map [#2424]
- Changed requirements to python-dateutil [#2426]
- Clarify coordinate system definitions [#2429]
- Improve Map Peek when using draw_grid [#2442]
- Add HCC --> HGS test [#2443]
- Testing the transformation linking SunPy and Astropy against published values [#2454]
- Fixed title bug in sunpy.timeseries.rhessi [#2477]
- Allow LineAnimator to accept a varying x-axis [#2491]
- Indexing Bug Fix to LineAnimator [#2560]
- Output sphinx warnings to stdout [#2553]
- Docstring improvement for LineAnimator [#2514]
- move the egg_info builds to circleci [#2512]
- Added tests for TraceMap [#2504]
- Fix HGS frame constructor and HPC ``calculate_distance`` with SkyCoord constructor. [#2463]
- removed `wavelnth` keyword in meta desc of Maps to avoid using non standard FITS keyword like `nan` [#2456]
- The documentation build now uses the Sphinx configuration from sphinx-astropy rather than from astropy-helpers.[#2494]
- Migrate to hypothesis.strategies.datetimes [#2368]
- Prevent a deprecation warning due to truth values of Quantity [#2358]
- Print a warning when heliographic longitude is set to it's default value of 0 [#2480]
- parse_time now parses numpy.datetime64 correctly. [#2572]

0.8.5
=====

Bug Fixes
---------

- Removed AstropyDeprecationWarning from sunpy.coordinates.representation [#2476]
- Fix for NorthOffsetFrame under Astropy 3.0 [#2486]
- Fix lightcurve tests under numpy dev [#2505]
- Updated depecration link of radiospectra [#2481]
- Fixed Padding values in some of the documentation pages [#2497]
- Move documentation build to circleci [#2509]
- Fix Issue #2470 hgs_to_hcc(heliogcoord, heliocframe) [#2502]
- Fixing CompositeMap object so that it respects masked maps [#2492]

0.8.4
=====

Bug Fixes
---------

- Improve detection of ``SkyCoord`` frame instantiation when distance is
  `1*u.one`. This fixes a plotting bug with ``WCSAxes`` in Astropy 3.0 [#2465]
- removed `wavelnth` keyword in meta desc of Maps to avoid using non standard FITS keyword like `nan` [#2427]
- Change the default units for HPC distance from `u.km` to `None`. [#2465]

0.8.3
=====

Bug Fixes
---------

- `~sunpy.net.dataretriever.clients.XRSClient` now reports time ranges of files correctly. [#2364]
- Make parse_time work with datetime64s and pandas series [#2370]
- CompositeMap axes scaling now uses map spatial units [#2310]
- Moved license file to root of repository and updated README file [#2326]
- Fix docstring formatting for net.vso.attrs [#2309]]
- Fix coloring of ticks under matplotlib 2.0 default style [#2320]
- Always index arrays with tuples in `ImageAnimator` [#2320]
- Added links to possible attrs for FIDO in guide [#2317] [#2289]
- Updated GitHub Readme [#2281] [#2283]
- Fix matplotlib / pandas 0.21 bug in examples [#2336]
- Fixes the off limb enhancement example [#2329]
- Changes to masking hot pixels and picking bright pixels examples [#2325] [#2319]
- Travis CI fix for numpy-dev build [#2340]
- Updated masking brightest pixel example [#2338]
- Changed TRAVIS cronjobs [#2338]
- Support array values for `obstime` for coordinates and transformations [#2342] [#2346]
- Updated Gallery off limb enhance example [#2337]
- Documentation fixes for VSO [#2354] [#2353]
- All tests within the documentation have been fixed [#2343]
- Change to using pytest-remotedata for our online tests [#2345]
- Fixed upstream astropy/numpy documentation issues [#2359]
- Documentation for Map improved [#2361]
- Fix the output units of pixel_to_world [#2362]
- Documentation for Database improved [#2355]
- Added test for mapsave [#2365]
- Documentation for Sun improved [#2369]

0.8.2
=====

Bug Fixes
---------

- Shows a warning if observation time is missing [#2293]
- Updates MapCube to access the correct properties of the namedtuple SpatialPair [#2297]

0.8.1
======

Bug fixes
---------

- Fixed TimeSeries test failures due to missing test files [#2273]
- Refactored a GOES test to avoid a Py3.6 issue [#2276]

0.8.0
======

New Features
------------

-  Solar differential rotation for maps and submaps included.
-  Solar rotation calculation and mapcube derotation now use sunpy coordinates.
-  Sample data now downloads automatically on import if not available
   and is now pluggable so can be used by affiliated packages. Shortcut
   names have been normalized and all LIGHTCURVE shortcuts have changed
   to TIMESERIES.
-  Calculation of points on an arc of a great circle connecting two
   points on the Sun.
-  Removed ``extract_time`` function from ``sunpy.time`` and also tests
   related to the function from ``sunpy.time.tests``
-  User can now pass a custom time format as an argument inside
   ``sunpy.database.add_from_dir()`` in case the ``date-obs`` metadata
   cannot be read automatically from the files.
-  Add time format used by some SDO HMI FITS keywords
-  Now the ``sunpy.database.tables.display_entries()`` prints an astropy
   table.
-  Additional methods added inside the ``sunpy.database`` class to make
   it easier to display the database contents.
-  Remove unused ``sunpy.visualization.plotting`` module
-  Port the pyana wrapper to Python 3
-  ``Map.peek(basic_plot-True)`` no longer issues warnings
-  Remove the ``sunpy.map.nddata_compat`` module, this makes
   ``Map.data`` and ``Map.meta`` read only.
-  Add a ``NorthOffsetFrame`` class for generating HGS-like coordinate
   systems with a shifted north pole.
-  Remove deprecated ``VSOClient.show`` method.
-  Deprecate ``sunpy.wcs``: ``sunpy.coordinates`` and ``sunpy.map`` now
   provide all that functionality in a more robust manner.
-  Added hdu index in ``sunpy.database.tables.DatabaseEntry`` as a
   column in the table.
-  Removed ``HelioviewerClient`` from the ``sunpy.net`` namespace. It
   should now be imported with
   ``from sunpy.net.helioviewer import HelioviewerClient``.
-  Removed compatibility with standalone ``wcsaxes`` and instead depend
   on the version in astropy 1.3. SunPy now therefore depends on
   astropy>-1.3.
-  Update to ``TimeRange.__repr__``; now includes the qualified name and
   ``id`` of the object.
-  A new ``sunpy.visualization.imageanimator.LineAnimator`` class has
   been added to animate 1D data. This has resulted in API change for
   the ``sunpy.visualization.imageanimator.ImageAnimator`` class. The
   updateimage method has been renamed to update\_plot.
-  Drop support for Python 3.4.
-  SunPy now requires WCSAxes and Map.draw\_grid only works with
   WCSAxes.
-  ``Helioprojective`` and ``HelioCentric`` frames now have an
   ``observer`` attribute which itself is a coordinate object
   (``SkyCoord``) instead of ``B0``, ``L0`` and ``D0`` to describe the
   position of the observer.
-  ``GenericMap.draw_grid`` now uses ``WCSAxes``, it will only work on a
   ``WCSAxes`` plot, this may be less performant than the previous
   implementation.
-  ``GenericMap.world_to_pixel`` and ``GenericMap.pixel_to_world`` now
   accept and return ``SkyCoord`` objects only.
-  ``GenericMap`` has a new property ``observer_coordinate`` which
   returns a ``SkyCoord`` describing the position of the observer.
-  ``GenericMap.submap`` now takes arguments of the form ``bottom_left``
   and ``top_right`` rather than ``range_a`` and ``range_b``. This
   change enables submap to properly handle rotated maps and take input
   in the form of ``SkyCoord`` objects.
-  When referring to physical coordinates ``Pair.x`` has been replaced
   with ``SpatialPair.axis1``. This means values returned by
   ``GenericMap`` now differentiate between physical and pixel
   coordinates.
-  The physical radius of the Sun (length units) is now passed from Map
   into the coordinate frame so a consistent value is used when
   calculating distance to the solar surface in the
   ``HelioprojectiveFrame`` coordinate frame.
-  A new ``sunpy.visualization.imageanimator.ImageAnimatorWCS`` class
   has been added to animate N-Dimensional data with the associated WCS
   object.
-  Moved Docs to docs/ to follow the astropy style
-  Added SunPy specific warnings under util.
-  SunPy coordinate frames can now be transformed to and from Astropy
   coordinate frames
-  The time attribute for SunPy coordinate frames has been renamed from
   ``dateobs`` to ``obstime``
-  Ephemeris calculations with higher accuracy are now available under
   ``sunpy.coordinates.ephemeris``
-  Add support for SunPy coordinates to specify observer as a string of
   a major solar-system body, with the default being Earth. To make
   transformations using an observer specified as a string, ``obstime``
   must be set.
-  Added VSO query result block level caching in the database module.
   This prevents re-downloading of files which have already been
   downloaded. Especially helpful in case of overlapping queries.
-  Change the default representation for the Heliographic Carrington
   frame so Longitude follows the convention of going from 0-360
   degrees.
-  All Clients that are able to search and download data now have a
   uniform API that is `search` and `fetch`. The older functions are
   still there but are deprecated for 0.8.

Bug fixes
---------

-  Add tests for RHESSI instrument
-  Maps from Helioviewer JPEG2000 files now have correct image scaling.
-  Get and set methods for composite maps now use Map plot\_settings.
-  Simplified map names when plotting.
-  Fix bug in ``wcs.convert_data_to_pixel`` where crpix[1] was used for
   both axes.
-  Fix some leftover instances of ``GenericMap.units``
-  Fixed bugs in ``sun`` equations
-  ``sunpy.io.fits.read`` will now return any parse-able HDUs even if
   some raise an error.
-  ``VSOClient`` no longer prints a lot of XML junk if the query fails.
-  Fix Map parsing of some header values to allow valid float strings
   like 'nan' and 'inf'.
-  Fix Map parsing of some header values to allow valid float strings
   like 'nan' and 'inf'.

0.7.8
=====

-  The SunPy data directory "~/sunpy" is no longer created until it is
   used (issue #2018)
-  Change the default representation for the Heliographic Carrington
   frame so Longitude follows the convention of going from 0-360
   degrees.
-  Fix for surface gravity unit.
-  Support for Pandas 0.20.1

0.7.7
=====

-  Fix errors with Numpy 1.12

0.7.6
=====

-  Add Astropy 1.3 Support

0.7.5
=====

-  Fix test faliure (mapbase) with 1.7.4
-  Restrict supported Astropy version to 1.0<astropy<1.3
-  Add Figure test env to SunPy repo.

0.7.4
=====

-  Remove Map always forcing warnings on.
-  ``Map.center`` now uses ``Map.wcs`` to correctly handle rotation.
-  Fix link in coordinates documentation.
-  Update helioviewer URL to HTTPS (fixes access to Helioviewer).
-  Fix processing of TRACE and YOHKOH measurement properties.
-  Remove warnings when using ``Map.peek(basic_plot-True)``
-  Update docstrings for HPC and HCC frames.

0.7.3
=====

-  Fix ConfigParser for Python 3.5.2 - This allows SunPy to run under
   Python 3.5.2
-  Fix incorrect ordering of keys in ``MapMeta``
-  Add ``sunpy.util.scraper`` to the API documentation.

0.7.2
=====

-  Fixed bugs in ``sun`` equations

0.7.1
=====

-  Fix bug in ``wcs.convert_data_to_pixel`` where crpix[1] was used for
   both axes.
-  Fix some leftover instances of ``GenericMap.units``
-  Fixed bugs in ``sun`` equations
-  Now the ``sunpy.database.tables.display_entries()`` prints an astropy
   table.
-  Additional methods added inside the ``sunpy.database`` class to make
   it easier to display the database contents.
-  ``sunpy.io.fits.read`` will now return any parse-able HDUs even if
   some raise an error.
-  ``VSOClient`` no longer prints a lot of XML junk if the query fails.
-  Remove unused ``sunpy.visualization.plotting`` module
-  ``Map.peek(basic_plot-True)`` no longer issues warnings
-  Remove the ``sunpy.map.nddata_compat`` module, this makes
   ``Map.data`` and ``Map.meta`` read only.
-  Add a ``NorthOffsetFrame`` class for generating HGS-like coordinate
   systems with a shifted north pole.
-  Remove deprecated ``VSOClient.show`` method.
-  Deprecate ``sunpy.wcs``: ``sunpy.coordinates`` and ``sunpy.map`` now
   provide all that functionality in a more robust manner.
-  Added hdu index in ``sunpy.database.tables.DatabaseEntry`` as a
   column in the table.
-  Removed ``HelioviewerClient`` from the ``sunpy.net`` namespace. It
   should now be imported with
   ``from sunpy.net.helioviewer import HelioviewerClient``.
-  Removed compatibility with standalone ``wcsaxes`` and instead depend
   on the version in astropy 1.3. SunPy now therefore depends on
   astropy>-1.3.
-  Update to ``TimeRange.__repr__``; now includes the qualified name and
   ``id`` of the object.
-  Change the default representation for the Heliographic Carrington
   frame so Longitude follows the convention of going from 0-360
   degrees.
-  Fix Map parsing of some header values to allow valid float strings
   like 'nan' and 'inf'.

0.7.0
=====

-  Fixed test failures with numpy developer version.[#1808]
-  Added ``timeout`` parameter in ``sunpy.data.download_sample_data()``
-  Fixed ``aiaprep`` to return properly sized map.
-  Deprecation warnings fixed when using image coalignment.
-  Sunpy is now Python 3.x compatible (3.4 and 3.5).
-  Added a unit check and warnings for map metadata.
-  Added IRIS SJI color maps.
-  Updated ``show_colormaps()`` with new string filter to show a subset
   of color maps.
-  Fixed MapCube animations by working around a bug in Astropy's
   ImageNormalize
-  Remove ``vso.QueryResponse.num_records()`` in favour of ``len(qr)``
-  Add a ``draw_rectangle`` helper to ``GenericMap`` which can plot
   rectangles in the native coordinate system of the map.
-  Added the ability to shift maps to correct for incorrect map
   location, for example.
-  Bug fix for RHESSI summary light curve values.
-  Mapcube solar derotation and coalignment now pass keywords to the
   routine used to shift the images, scipy.ndimage.interpolation.shift.
-  Add automatic registration of ``GenericMap`` subclasses with the
   factory as long as they define an ``is_datasource_for`` method.
-  Added functions ``flareclass_to_flux`` and ``flux_to_flareclass``
   which convert between GOES flux to GOES class numbers (e.g. X12,
   M3.4).
-  Removed old ``sunpy.util.goes_flare_class()``
-  Bug fix for RHESSI summary light curve values.
-  The ``MapCube.as_array`` function now returns a masked numpy array if
   at least one of the input maps in the MapCube has a mask.
-  Map superpixel method now respects maps that have masks.
-  Map superpixel method now accepts numpy functions as an argument, or
   any user-defined function.
-  Map superpixel method no longer has the restriction that the number
   of original pixels in the x (or y) side of the superpixel exactly
   divides the number of original pixels in the x (or y) side of the
   original map data.
-  ``sunpy.physics.transforms`` has been deprecated and the code moved
   into ``sunpy.physics``.
-  Add the ``sunpy.coordinates`` module, this adds the core physical
   solar coordinates frame within the astropy coordinates framework.
-  Added ability of maps to draw contours on top of themselves
   (``draw_contours``)
-  Added concatenate functionality to lightcurve base class.
-  Fix Map to allow astropy.io.fits Header objects as valid input for
   meta arguments.
-  Added an examples gallery using ``sphinx-gallery``.
-  API clean up to constants. Removed constant() function which is now
   replaced by get().
-  Prevent helioviewer tests from checking access to the API endpoint
   when running tests offline.
-  ``GenericMap.units`` is renamed to ``GenericMap.spatial_units`` to
   avoid confusion with ``NDData.unit``.
-  ``GenericMap`` now has a ``coordinate_frame`` property which returns
   an ``astropy.coordinates`` frame with all the meta data from the map
   populated.
-  ``GenericMap`` now has a ``_mpl_axes`` method which allows it to be
   specified as a projection to ``matplotlib`` methods and will return a
   ``WCSAxes`` object with ``WCS`` projection.

0.6.5
=====

-  The draw\_grid keyword of the peek method of Map now accepts booleans
   or astropy quantities.
-  Fix bug in ``wcs.convert_data_to_pixel`` where crpix[1] was used for
   both axes.
-  Fixed bugs in ``sun`` equations

0.6.4
=====

-  Bug fix for rhessi summary lightcurve values.
-  Fix docstring for ``pixel_to_data`` and ``data_to_pixel``.
-  Fix the URL for the Helioviewer API. (This fixes Helioviewer.)
-  Fix the way ``reshape_image_to_4d_superpixel`` checks the dimension
   of the new image.
-  Fix Map to allow astropy.io.fits Header objects as valid input for
   meta arguments.
-  Prevent helioviewer tests from checking access to API when running
   tests in offline mode.

0.6.3
=====

-  Change setup.py extras to install suds-jurko not suds.

0.6.2
=====

-  Changed start of GOES 2 operational time range back to 1980-01-04 so
   data from 1980 can be read into GOESLightCurve object
-  Fix bug with numpy 1.10
-  update astropy\_helpers
-  Added new sample data

0.6.1
=====

-  Fixed MapCube animations by working around a bug in Astropy's
   ImageNormalize
-  Small fix to RTD builds for Affiliated packages
-  SunPy can now be installed without having to install Astropy first.
-  MapCubes processed with ``coalignment.apply_shifts`` now have correct
   metadata.
-  Multiple fixes for WCS transformations, especially with solar-x,
   solar-y CTYPE headers.

0.6.0
=====

-  Enforced the use of Astropy Quantities through out most of SunPy.
-  Dropped Support for Python 2.6.
-  Remove old style string formatting and other 2.6 compatibility lines.
-  Added vso like querying feature to JSOC Client.
-  Refactor the JSOC client so that it follows the .query() .get()
   interface of VSOClient and UnifedDownloader.
-  Provide ``__str__`` and ``__repr__`` methods on vso ``QueryResponse``
   deprecate ``.show()``.
-  Downloaded files now keep file extensions rather than replacing all
   periods with underscores.
-  Update to TimeRange API, removed t1 and t0, start and end are now
   read-only attributes.
-  Added ability to download level3 data for lyra Light Curve along with
   corresponding tests.
-  Added support for gzipped FITS files.
-  Add STEREO HI Map subclass and color maps.
-  Map.rotate() no longer crops any image data.
-  For accuracy, default Map.rotate() transformation is set to
   bi-quartic.
-  ``sunpy.image.transform.affine_transform`` now casts integer data to
   float64 and sets NaN values to 0 for all transformations except
   scikit-image rotation with order <- 3.
-  CD matrix now updated, if present, when Map pixel size is changed.
-  Removed now-redundant method for rotating IRIS maps since the
   functionality exists in Map.rotate()
-  Provide ``__str__`` and ``__repr__`` methods on vso ``QueryResponse``
   deprecate ``.show()``
-  SunPy colormaps are now registered with matplotlib on import of
   ``sunpy.cm``
-  ``sunpy.cm.get_cmap`` no longer defaults to 'sdoaia94'
-  Added database url config setting to be setup by default as a sqlite
   database in the sunpy working directory
-  Added a few tests for the sunpy.roi module
-  Added capability for figure-based tests
-  Removed now-redundant method for rotating IRIS maps since the
   functionality exists in Map.rotate().
-  SunPy colormaps are now registered with matplotlib on import of
   ``sunpy.cm``.
-  ``sunpy.cm.get_cmap`` no longer defaults to 'sdoaia94'.
-  Added database url config setting to be setup by default as a sqlite
   database in the sunpy working directory.
-  Added a few tests for the sunpy.roi module.
-  Refactored mapcube co-alignment functionality.
-  Removed sample data from distribution and added ability to download
   sample files
-  Changed start of GOES 2 operational time range back to 1980-01-04 so
   data from 1980 can be read into GOESLightCurve object
-  Require JSOC request data calls have an email address attached.
-  Calculation of the solar rotation of a point on the Sun as seen from
   Earth, and its application to the de-rotation of mapcubes.
-  Downloaded files now keep file extensions rather than replacing all
   periods with underscores
-  Fixed the downloading of files with duplicate names in sunpy.database
-  Removed sample data from distribution and added ability to download
   sample files.
-  Added the calculation of the solar rotation of a point on the Sun as
   seen from Earth, and its application to the de-rotation of mapcubes.
-  Changed default for GOESLightCurve.create() so that it gets the data
   from the most recent existing GOES fits file.
-  Map plot functionality now uses the mask property if it is present,
   allowing the plotting of masked map data
-  Map Expects Quantities and returns quantities for most parameters.
-  Map now used Astropy.wcs for world <-> pixel conversions.
-  map.world\_to\_pixel now has a similar API to map.pixel\_to\_world.
-  map.shape has been replaced with map.dimensions, which is ordered x
   first.
-  map.rsun\_arcseconds is now map.rsun\_obs as it returns a quantity.
-  Map properties are now named tuples rather than dictionaries.
-  Improvement for Map plots, standardization and improved color tables,
   better access to plot variables through new plot\_settings variable.
-  Huge improvements in Instrument Map doc strings. Now contain
   instrument descriptions as well as reference links for more info.
-  net.jsoc can query data series with time sampling by a Sample
   attribute implemented in vso.
-  MapCube.plot and MapCube.peek now support a user defined
   plot\_function argument for customising the animation.
-  Added new sample data file, an AIA cutout file.
-  Moved documentation build directory to doc/build

0.5.5
=====

-  Changed default for GOESLightCurve.create() so that it gets the data
   from the most recent existing GOES fits file.
-  Improvements to the Map documentation.
-  Typo fixes in sunpy.wcs documentation.

0.5.4
=====

-  ``sunpy.image.transform.affine_transform`` now casts integer data to
   float64 and sets NaN values to 0 for all transformations except
   scikit-image rotation with order <- 3.
-  Updated SWPC/NOAA links due to their new website.
-  Exposed the raw AIA color tables in ``sunpy.cm.color_tables``.
-  Fixes ``map`` compatibility with Astropy 1.0.x.

0.5.3
=====

-  Goes peek() plot now works with matplotlib 1.4.x
-  The ANA file reading C extensions will no longer compile under
   windows. Windows was not a supported platform for these C extensions
   previously.

0.5.2
=====

-  If no CROTA keyword is specified in Map meta data, it will now
   default to 0 as specified by the FITS WCS standard.
-  Map now correctly parses and converts the CD matrix, as long as CDELT
   is specified as well. (Fixes SWAP files)
-  Fix of HELIO webservice URLs
-  MapCube.plot() is now fixed and returns a
   matplotlib.animation.FuncAnimation object.

0.5.1
=====

-  MAJOR FIX: map.rotate() now works correctly for all submaps and off
   center rotations.
-  HELIO URL updated, querys should now work as expected.
-  All tabs removed from the code base.
-  All tests now use tempfile rather than creating files in the current
   directory.
-  Documentation builds under newer sphinx versions.
-  ANA and JP2 tests are skipped if dependancies are missing.
-  ANA tests are skipped on windows.

0.5.0
=====

-  Added additional functionality to the GOES module i.e. the ability to
   calculate GOES temperature and emission measure from GOES fluxes.
-  changed \_maps attribute in MapCube to a non-hidden type
-  Added Nobeyama Radioheliograph data support to Lightcurve object.
-  Fixed some tests on map method to support Windows
-  Added a window/split method to time range
-  Updates to spectrogram documentation
-  Added method Database.add\_from\_hek\_query\_result to HEK database
-  Added method Database.download\_from\_vso\_query\_result
-  GOES Lightcurve now makes use of a new source of GOES data, provides
   metadata, and data back to 1981.
-  Removed sqlalchemy as a requirement for SunPy
-  Added support for NOAA solar cycle prediction in lightcurves
-  Some basic tests for GenericLightCurve on types of expected input.
-  Fix algorithm in sunpy.sun.equation\_of\_center
-  Added Docstrings to LightCurve methods.
-  Added tests for classes in sunpy.map.sources. Note that some classes
   (TRACE, RHESSI) were left out because SunPy is not able to read their
   FITS files.
-  Added functions that implement image coalignment with support for
   MapCubes.
-  Cleaned up the sunpy namespace, removed .units, /ssw and .sphinx.
   Also moved .coords .physics.transforms.
-  Added contains functionality to TimeRange module
-  Added t-'now' to parse\_time to privide utcnow datetime.
-  Fixed time dependant functions (.sun) to default to t-'now'
-  Fixed solar\_semidiameter\_angular\_size
-  Improved line quality and performances issues with map.draw\_grid()
-  Remove deprecated make\_map command.

0.4.2
=====

-  Fixes to the operational range of GOES satellites
-  Fix the URL for HELIO queries.

0.4.1
=====

-  Fix map.rotate() functionality
-  Change of source for GOES data.
-  Fix EIT test data and sunpy FITS saving
-  Some documentation fixes
-  fix file paths to use os.path.join for platform independance.

0.4.0
=====

-  **Major** documentation refactor. A far reaching re-write and
   restructure.
-  Add a SunPy Database to store and search local data.
-  Add beta support for querying the HELIO HEC
-  Add beta HEK to VSO query translation.
-  Add the ability to download the GOES event list.
-  Add support for downloading and querying the LYTAF database.
-  Add support for ANA data.
-  Updated sun.constants to use astropy.constants objects which include
   units, source, and error instide. For more info check out
   http://docs.astropy.org/en/latest/constants/index.html
-  Add some beta support for IRIS data products
-  Add a new MapCubeAnimator class with interactive widgets which is
   returned by mapcube.peek().
-  The Glymur library is now used to read JPEG2000 files.
-  GOESLightCurve now supports all satellites.
-  Add support for VSO queries through proxies.
-  Fix apparent Right Ascension calulations.
-  LightCurve meta data member now an OrderedDict Instance

0.3.2
=====

-  Pass draw\_limb arguments to patches.Circle
-  Pass graw\_grid arguments to pyplot.plot()
-  Fix README code example
-  Fix Documentation links in potting guide
-  Update to new EVE data URL
-  Update LogicalLightcurve example in docs
-  Improved InteractiveVSOClient documentation
-  GOESLightCurve now fails politely if no data is avalible.

Known Bugs:

-  sunpy.util.unit\_conversion.to\_angstrom does not work if 'nm' is
   passed in.

0.3.1
=====

-  Bug Fix: Fix a regression in CompositeMap that made contor plots
   fail.
-  Bug Fix: Allow Map() to accept dict as metadata.
-  Bug Fix: Pass arguments from Map() to io.read\_file.

0.3.0
=====

-  Removal of Optional PIL dependancy
-  Parse\_time now looks through nested lists/tuples
-  Draw\_limb and draw\_grid are now implemented on MapCube and
   CompositeMap
-  Caculations for differential roation added
-  mapcube.plot() now runs a mpl animation with optional controls
-  A basic Region of Interest framework now exists under sunpy.roi
-  STEREO COR colour maps have been ported from solarsoft.
-  sunpy.time.timerange has a split() method that divides up a time
   range into n equal parts.
-  Added download progress bar
-  pyfits is depricated in favor of Astropy

spectra:

-  Plotting has been refactorted to use a consistent interface
-  spectra now no-longer inherits from numpy.ndarray instead has a .data
   attribute.

Map: \* map now no-longer inherits from numpy.ndarray instead has a
.data attribute. \* make\_map is deprecated in favor of Map which is a
new factory class \* sunpy.map.Map is now sunpy.map.GenericMap \*
mymap.header is now mymap.meta \* attributes of the map class are now
read only, changes have to be made through map.meta \* new MapMeta class
to replace MapHeader, MapMeta is not returned by sunpy.io \* The
groundwork for GenericMap inherting from astropy.NDData has been done,
there is now a NDDataStandin class to provide basic functionality.

io: \* top level file\_tools improved to be more flexible and support
multiple HDUs \* all functions in sunpy.io now assume mutliple HDUs,
even JP2 ones. \* there is now a way to override the automatic filetype
detection \* Automatic fits file detection improved \* extract\_waveunit
added to io.fits for detection of common ways of storing wavelength unit
in fits files.

-  A major re-work of all interal imports has resulted in a much cleaner
   namespace, i.e. sunpy.util.util is no longer used to import util.
-  Some SOHO and STEREO files were not reading properly due to a
   date\_obs parameter.
-  Sunpy will now read JP2 files without a comment parameter.
-  Memory leak in Crotate patched
-  Callisto: Max gap between files removed

0.2.0
=====

-  Completely re-written plotting routines for most of the core
   datatypes.
-  JPEG 2000 support as an input file type.
-  Improved documentation for much of the code base, including
   re-written installation instructions.
-  New lightcurve object

   -  LYRA support
   -  GOES/XRS support
   -  SDO/EVE support

-  New Spectrum and Spectrogram object (in development)

   -  Spectrogram plotting routines
   -  Callisto spectrum type and support
   -  STEREO/SWAVES support

-  Map Object

   -  Added support for LASCO, Yohkoh/XRT maps
   -  A new CompositeMap object for overlaying maps
   -  Resample method
   -  Superpixel method
   -  The addition of the rotate() method for 2D maps.
