Sunpy v1.0.9 (2020-03-27)
=========================

Bug Fixes
---------

- Fix a bug in `sunpy.net.jsoc.JSOCClient` where requesting data for export would not work if a non-time primekey was used. (`#3825 <https://github.com/sunpy/sunpy/pull/3825>`__)
- Add explicit support for dealing with download urls for files, under 'as-is' protocol in `sunpy.net.jsoc.JSOCClient.get_request`. (`#3838 <https://github.com/sunpy/sunpy/pull/3838>`__)
- Add support for passing paths of type `pathlib.Path` in `sunpy.net.jsoc.JSOCClient.fetch`. (`#3838 <https://github.com/sunpy/sunpy/pull/3838>`__)
- Fix failing of fetching of the indexed JSOCResponses using `Fido.fetch`. (`#3852 <https://github.com/sunpy/sunpy/pull/3852>`__)
- Prevented `GenericMap.plot` modifying in-place any items passed as ``imshow_kwargs``. (`#3867 <https://github.com/sunpy/sunpy/pull/3867>`__)
- Changed the format of DATE-OBS in `GenericMap.wcs` from iso to isot (ie. with a "T" between the date and time) to conform with the FITS standard. (`#3872 <https://github.com/sunpy/sunpy/pull/3872>`__)
- Fixed a minor error (up to ~10 arcseconds) in the calculation of the Sun's position angle (:func:`sunpy.coordinates.sun.P`). (`#3886 <https://github.com/sunpy/sunpy/pull/3886>`__)
- `~sunpy.net.hek.HEKClient` was returning HTML and not JSON. (`#3899 <https://github.com/sunpy/sunpy/pull/3899>`__)
- Updated to HTTPS for HEK. (`#3917 <https://github.com/sunpy/sunpy/pull/3917>`__)


Improved Documentation
----------------------

- Changed padding value of an example in the example gallery to fix the overlap of titles and x-label axes. (`#3835 <https://github.com/sunpy/sunpy/pull/3835>`__)
- Clarified some inputs to `sunpy.map.GenericMap.plot`. (`#3866 <https://github.com/sunpy/sunpy/pull/3866>`__)


Trivial/Internal Changes
------------------------

- Created a helper function for testing the equality/closeness of longitude angles (i.e., angles with wrapping). (`#3804 <https://github.com/sunpy/sunpy/pull/3804>`__)


Sunpy v1.0.8 (2020-02-13)
=========================

Features
--------

- Updated the gallery example titled 'Downloading and plotting an HMI magnetogram' to rotate the HMI magnetogram such that solar North is pointed up. (`#3573 <https://github.com/sunpy/sunpy/pull/3573>`__)


Bug Fixes
---------

- Fixed a bug where permission denied errors when downloading files are very verbose by adding an error message in `~sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch`. (`#3417 <https://github.com/sunpy/sunpy/pull/3417>`__)
- Added support for passing ``TimeSeriesMetaData`` object to ``timeseries_factory`` and associated validation tests. (`#3639 <https://github.com/sunpy/sunpy/pull/3639>`__)
- Work around incorrect Content-Disposition headers in some VSO downloads, which were leading to mangled filenames. (`#3740 <https://github.com/sunpy/sunpy/pull/3740>`__)
- Fixed a bug with the calculation of Carrington longitude as seen from Earth where it was using an old approach instead of the current approach (for example, the varying Sun-Earth distance is now taken into account).
  The old approach resulted in errors no greater than 7 arcseconds in Carrington longitude when using `~sunpy.coordinates.sun.L0` and `~sunpy.coordinates.frames.HeliographicCarrington`. (`#3772 <https://github.com/sunpy/sunpy/pull/3772>`__)


Improved Documentation
----------------------

- Fixed the plots with multiple subplots in the ``Map`` user guide to properly use `~astropy.visualization.wcsaxes` and to be appropriately sized. (`#3454 <https://github.com/sunpy/sunpy/pull/3454>`__)
- A new example gallery example "Plotting a difference image" has been added,
  which can be used for base difference or running difference images. (`#3627 <https://github.com/sunpy/sunpy/pull/3627>`__)
- Corrected misleading `~sunpy.timeseries.metadata.TimeSeriesMetaData` documentation about optional parameters. (`#3680 <https://github.com/sunpy/sunpy/pull/3680>`__)


Trivial/Internal Changes
------------------------

- Fixed the transformation test for `~sunpy.coordinates.metaframes.NorthOffsetFrame`, which would intermittently fail. (`#3775 <https://github.com/sunpy/sunpy/pull/3775>`__)


Sunpy v1.0.7 (2020-01-10)
=========================

Bug Fixes
---------

- Fixed bugs with some coordinate transformations when ``obstime`` is ``None`` on the destination frame but can be assumed to be the same as the ``obstime`` of the source frame. (`#3576 <https://github.com/sunpy/sunpy/pull/3576>`__)
- Updated `sunpy.map.mapsequence.MapSequence` so that calling ``_derotate()`` raises ``NotImplementedError``.
  Added associated tests. (`#3613 <https://github.com/sunpy/sunpy/pull/3613>`__)
- Fixed pandas plotting registration in `sunpy.timeseries`. (`#3633 <https://github.com/sunpy/sunpy/pull/3633>`__)


Improved Documentation
----------------------

- Clarified the meaning of some `GenericMap` observer properties. (`#3585 <https://github.com/sunpy/sunpy/pull/3585>`__)
- Added inherited members of `sunpy.map` classes to the docs. (`#3587 <https://github.com/sunpy/sunpy/pull/3587>`__)
- Fixed documentation of `sunpy.database.Database.search` by adding ``Returns`` docstring. (`#3593 <https://github.com/sunpy/sunpy/pull/3593>`__)
- Updated the docstring for the parameter ``sortby`` in `~sunpy.map.MapSequence` with the default value, valid value and how to disable sorting. (`#3601 <https://github.com/sunpy/sunpy/pull/3601>`__)
- Updated the tour guide to reflect that the time series is not random data. (`#3603 <https://github.com/sunpy/sunpy/pull/3603>`__)


Sunpy v1.0.6 (2019-11-20)
=========================

Bug Fixes
---------

- `~sunpy.coordinates.utils.GreatArc` now accounts for the start and end points of the arc having different observers. (`#3334 <https://github.com/sunpy/sunpy/pull/3334>`__)
- Single character wildcards and character ranges can now be passed as
  glob patterns to `~sunpy.map.Map`. (`#3408 <https://github.com/sunpy/sunpy/pull/3408>`__)
- The returned list of `~sunpy.map.Map` objects is now sorted by filename when
  passing a directory or glob pattern to `~sunpy.map.MapFactory`. (`#3408 <https://github.com/sunpy/sunpy/pull/3408>`__)
- Fixed a bug where clipping behavior had been enabled by default in the plotting normalizers for ``Map`` objects.  Clipping needs to be disabled to make use of the over/under/masked colors in the colormap. (`#3427 <https://github.com/sunpy/sunpy/pull/3427>`__)
- Fix a bug with observer based frames that prevented a coordinate with an array of obstimes being transformed to other frames. (`#3455 <https://github.com/sunpy/sunpy/pull/3455>`__)
- `sunpy.map.GenericMap` will no longer raise a warning if the posisition of the
  observer is not known for frames that don't need an observer, i.e. heliographic
  frames. (`#3462 <https://github.com/sunpy/sunpy/pull/3462>`__)
- Apply `os.path.expanduser` to `sunpy.map.MapFactory` input
  before passing to `glob.glob` (`#3477 <https://github.com/sunpy/sunpy/pull/3477>`__)
- Fix multiple instances of `sunpy.map.sources` assuming the type of FITS Header
  values. (`#3497 <https://github.com/sunpy/sunpy/pull/3497>`__)


Improved Documentation
----------------------

- Clarified the meaning of :attr:`GenericMap.dsun`. (`#3430 <https://github.com/sunpy/sunpy/pull/3430>`__)
- Updated the user guide for Map to use ``clip_interval``. (`#3450 <https://github.com/sunpy/sunpy/pull/3450>`__)
- Updated the Venus-transit gallery to use the VSO so that it has correct pointing information in the header. (`#3451 <https://github.com/sunpy/sunpy/pull/3451>`__)
- Fixed various issues with the gallery example of saving/loading coordinates using `asdf`. (`#3473 <https://github.com/sunpy/sunpy/pull/3473>`__)
- Added ``sunpy.__citation__`` with a BibTex entry for citing sunpy. (`#3478 <https://github.com/sunpy/sunpy/pull/3478>`__)
- Added an example showing how to display two maps and fade between them. (`#3488 <https://github.com/sunpy/sunpy/pull/3488>`__)


Trivial/Internal Changes
------------------------

- Copy the library `distro` into `sunpy/extern`: replaces the deprecated `platform/linux_distribution` (`#3396 <https://github.com/sunpy/sunpy/pull/3396>`__)
- Corrected spelling of 'plotting' in timeseries method (changed 'ploting' to 'plotting'). (`#3429 <https://github.com/sunpy/sunpy/pull/3429>`__)


Sunpy v1.0.5 (2019-10-22)
=========================

Bug Fixes
---------

- Fix incorrect files being included in the tarball, and docs missing from the
  tarball (`#3423 <https://github.com/sunpy/sunpy/pull/3423>`__)


Sunpy v1.0.4 (2019-10-22)
=========================

Bug Fixes
---------

- Fixed situations where 2D coordinates provided to `~sunpy.coordinates.frames.HeliographicStonyhurst` and `~sunpy.coordinates.frames.HeliographicCarrington` were not converted to 3D as intended.  Furthermore, the stored data will always be the post-conversion, 3D version. (`#3351 <https://github.com/sunpy/sunpy/pull/3351>`__)
- Fix off by one error in `sunpy.map.make_fitswcs_header` where when using the
  default ``reference_pixel=None`` keyword argument the pixel coordinate of the
  reference pixel was off by +1. (`#3356 <https://github.com/sunpy/sunpy/pull/3356>`__)
- Fixing the ordering of lon and lat inputs into `sunpy.map.make_fitswcs_header` (`#3371 <https://github.com/sunpy/sunpy/pull/3371>`__)
- Updated the URL for Fermi spacecraft-pointing files to use an HTTPS connection to HEASARC. (`#3381 <https://github.com/sunpy/sunpy/pull/3381>`__)


Improved Documentation
----------------------

- Improved the contributing guide by updating commands and highlighting text. (`#3394 <https://github.com/sunpy/sunpy/pull/3394>`__)
- Removing `.fits` from the end of path kwargs in `sunpy.net.FIDO.fetch` docs to change output file extension from `{file}.fits.fits` to `{file}.fits`. (`#3399 <https://github.com/sunpy/sunpy/pull/3399>`__)


Sunpy v1.0.3 (2019-08-29)
=========================

Features
--------

- Add ability to disable progressbars when dowloading files using `sunpy.net.helioviewer.py` and edited docstrings to mention this feature. (`#3280 <https://github.com/sunpy/sunpy/pull/3280>`__)


Bug Fixes
---------

- Fixed the handling of coordinates with velocity information when transforming between Astropy frames and SunPy frames. (`#3247 <https://github.com/sunpy/sunpy/pull/3247>`__)
- Fixed all coordinate transformations to properly handle a change in observation time. (`#3247 <https://github.com/sunpy/sunpy/pull/3247>`__)
- Fixed `~sunpy.physics.solar_rotation.calculate_solar_rotate_shift` so that it does not calculate a shift between the reference layer and itself, which would sometimes incorrectly result in a shift of a pixel due to numerical precision. (`#3255 <https://github.com/sunpy/sunpy/pull/3255>`__)
- Stop crash when ``LineAnimator`` ``axes_ranges`` entry given as ``1D`` array when data is ``>1D``, i.e. as an independent axis. (`#3283 <https://github.com/sunpy/sunpy/pull/3283>`__)
- Fixed a bug where the transformation from `~sunpy.coordinates.frames.Helioprojective` to `~sunpy.coordinates.frames.Heliocentric` used the Sun-observer distance from the wrong frame when shifting the origin, and thus might not give the correct answer if the observer was not the same for the two frames. (`#3291 <https://github.com/sunpy/sunpy/pull/3291>`__)
- Fixed a bug with the transformations between `~sunpy.coordinates.frames.Heliocentric` and `~sunpy.coordinates.frames.HeliographicStonyhurst` when the frame observation time was not the same as the observer observation time.  The most common way to encounter this bug was when transforming from `~sunpy.coordinates.frames.Helioprojective` to any non-observer-based frame while also changing the observation time. (`#3291 <https://github.com/sunpy/sunpy/pull/3291>`__)
- Fixed a `sunpy.coordinates` bug where a frame using the default observer of Earth could have its observer overwritten during a transformation. (`#3291 <https://github.com/sunpy/sunpy/pull/3291>`__)
- VSO client `fetch` should not download when `wait` keyword argument is specified. (`#3298 <https://github.com/sunpy/sunpy/pull/3298>`__)
- Fixed a bug with `~sunpy.coordinates.wcs_utils.solar_frame_to_wcs_mapping` that assumed that the supplied frame was a SunPy frame. (`#3305 <https://github.com/sunpy/sunpy/pull/3305>`__)
- Fixed bugs with `~sunpy.coordinates.wcs_utils.solar_frame_to_wcs_mapping` if the input frame does not include an observation time or an observer. (`#3305 <https://github.com/sunpy/sunpy/pull/3305>`__)


Improved Documentation
----------------------

- Added more details to docstrings in `sunpy.coordinates.frames`. (`#3262 <https://github.com/sunpy/sunpy/pull/3262>`__)
- Added a link to package maintainer list in the API Stability page. (`#3281 <https://github.com/sunpy/sunpy/pull/3281>`__)


Trivial/Internal Changes
------------------------

- Allow running our sphinx-gallery examples as Jupyter notebooks via Binder (`#3256 <https://github.com/sunpy/sunpy/pull/3256>`__)


Sunpy v1.0.2 (2019-06-26)
=========================

Bug Fixes
---------

- `sunpy.map.sources.AIAMap` and `sunpy.map.sources.HMIMap` will no longer assume
  the existance of certain header keys. (`#3217 <https://github.com/sunpy/sunpy/pull/3217>`__)
- `sunpy.map.make_fitswcs_header` now supports specifying the map projection
  rather than defaulting to ``TAN``. (`#3218 <https://github.com/sunpy/sunpy/pull/3218>`__)
- Fix the behaviour of
  `sunpy.coordinates.frames.Helioprojective.calculate_distance` if the
  representation isn't Spherical. (`#3219 <https://github.com/sunpy/sunpy/pull/3219>`__)
- Fixed a bug where the longitude of a coordinate would not wrap at the expected angle following a frame transformation. (`#3223 <https://github.com/sunpy/sunpy/pull/3223>`__)
- Fixed a bug where passing a time or time interval to the differential rotation function threw an error because the new observer was not in HGS. (`#3225 <https://github.com/sunpy/sunpy/pull/3225>`__)
- Fixed bug where `~sunpy.coordinates.ephemeris.get_horizons_coord` was unable to accept `~astropy.time.Time` arrays as input. (`#3227 <https://github.com/sunpy/sunpy/pull/3227>`__)
- Fix the ticks on the default heliographic grid overlay so they are not white
  (and normally invisible) by default. (`#3235 <https://github.com/sunpy/sunpy/pull/3235>`__)
- Fixed a bug with `sunpy.net.hek.HEKClient` when the results returned were a mixed dataset. (`#3240 <https://github.com/sunpy/sunpy/pull/3240>`__)
- Fix `sunpy.physics.differential_rotation.differential_rotate` to rotate in the
  correct direction and to account for the rotation of the heliographic
  coordinate frame with time. (`#3245 <https://github.com/sunpy/sunpy/pull/3245>`__)
- Fixed a bug with the handling of changing observation times for transformations between `~astropy.coordinates.HCRS` and `~sunpy.coordinates.frames.HeliographicStonyhurst`, which also indirectly affected other transformations when changing observation times. (`#3246 <https://github.com/sunpy/sunpy/pull/3246>`__)


Improved Documentation
----------------------

- Cleaned and expanded upon the docstrings for each Fido Client. (`#3220 <https://github.com/sunpy/sunpy/pull/3220>`__)
- Added clarifying hyperlinks to the gallery example `getting_lasco_observer_location` to link to `astroquery` docs page. (`#3228 <https://github.com/sunpy/sunpy/pull/3228>`__)


Sunpy v1.0.1 (2019-06-07)
=========================

Bug Fixes
---------

- Fixed accuracy issues with the calculations of Carrington longitude (`~sunpy.coordinates.sun.L0`) and Carrington rotation number (`~sunpy.coordinates.sun.carrington_rotation_number`). (`#3178 <https://github.com/sunpy/sunpy/pull/3178>`__)
- Updated `sunpy.map.header_helper.make_fitswcs_header` to be more strict on the inputs it accepts. (`#3183 <https://github.com/sunpy/sunpy/pull/3183>`__)
- Fix the calculation of ``rsun_ref`` in `~sunpy.map.make_fitswcs_header` and and
  ensure that the default reference pixel is indexed from 1. (`#3184 <https://github.com/sunpy/sunpy/pull/3184>`__)
- Fixed the missing transformation between two `~sunpy.coordinates.HeliographicCarrington` frames with different observation times. (`#3186 <https://github.com/sunpy/sunpy/pull/3186>`__)


Improved Documentation
----------------------

- Clean up the docstring for `sunpy.physics.differential_rotation.solar_rotate_coordinate` to make the example clearer. (`#2708 <https://github.com/sunpy/sunpy/pull/2708>`__)
- Added new gallery examples and cleaned up various gallery examples. (`#3181 <https://github.com/sunpy/sunpy/pull/3181>`__)


Sunpy 1.0.0 (2019-06-01)
========================

Backwards Incompatible Changes
------------------------------

- Move the matplotlib animators from ``sunpy.visualisation.imageanimator`` and
  ``sunpy.visualization.mapcubeanimator`` to `sunpy.visualization.animator`. (`#2515 <https://github.com/sunpy/sunpy/pull/2515>`__)
- Make `sunpy.time.parse_time` return `astropy.time.Time` instead of `datetime.datetime`. (`#2611 <https://github.com/sunpy/sunpy/pull/2611>`__)
- The properties and methods of `sunpy.time.TimeRange` returns `astropy.time.Time` and `astropy.time.TimeDelta` instead of `datetime.datetime` and `datetime.timedelta` respectively. (`#2638 <https://github.com/sunpy/sunpy/pull/2638>`__)
- The `sunpy.instr.goes` module now accepts and returns
  `sunpy.timeseries.XRSTimeSeries` objects only. (`#2666 <https://github.com/sunpy/sunpy/pull/2666>`__)
- ``obstime`` keyword param of ``sunpy.instr.goes._goes_lx`` takes a non-scalar `astropy.time.Time` object instead of `numpy.ndarray`. The precision of times contained in `sunpy.timeseries` has been increased to 9 from 6. (`#2676 <https://github.com/sunpy/sunpy/pull/2676>`__)
- Removed ``sunpy.net.jsoc.attrs.Time`` because it served the same purpose as `sunpy.net.attrs.Time` after the switch to `astropy.time.Time`. (`#2694 <https://github.com/sunpy/sunpy/pull/2694>`__)
- Remove unused ``**kwargs`` within TimeSeries functions. (`#2717 <https://github.com/sunpy/sunpy/pull/2717>`__)
- Rotation matrices inside map objects were previously stored as numpy matrices, but are now
  stored as numpy arrays, as numpy will eventually remove their matrix datatype. See
  https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html for more information. (`#2719 <https://github.com/sunpy/sunpy/pull/2719>`__)
- The `sunpy.cm.show_colormaps` function now accepts the keyword 'search' instead of 'filter'. (`#2731 <https://github.com/sunpy/sunpy/pull/2731>`__)
- The keyword arguments to all client ``.fetch`` methods have been changed to
  support the new parfive downloader and to ensure consisteny across all Fido
  clients. (`#2797 <https://github.com/sunpy/sunpy/pull/2797>`__)
- The Helioviewer client has been switched to using the newer Helioviewer API.
  This has meant that we have changed some of the keywords that were passed into client's methods.
  We have enforced that several keywords (observatory,instrument,detector,measurement) need to be defined otherwise the functions cannot return any data. (`#2801 <https://github.com/sunpy/sunpy/pull/2801>`__)
- Maps no longer assume that the pixel units are arcseconds if the units aren't
  explicitly set. In addition to this if critical metadata is missing from when
  creating a map, the map will fail to initialize and will raise an error. (`#2847 <https://github.com/sunpy/sunpy/pull/2847>`__)
- axis_ranges kwarg of `sunpy.visualization.animator.base.ArrayAnimator`, `sunpy.visualization.animator.image.ImageAnimator` and `sunpy.visualization.animator.line.LineAnimator` now must be entered as None, [min, max] or pixel edges of each array element. Previously, pixel centers were expected.  This change removes ambiguity in interpretation and ensures the extent of the plot can always be accurately derived. (`#2867 <https://github.com/sunpy/sunpy/pull/2867>`__)
- All keywords have been added (with defaults) to each `~sunpy.net.helioviewer.HelioviewerClient` function.
  This means that there will be some changes to the style of the PNG screenshot that is returned.
  Returns for the JPEG 2000 and the other functions should be the same but not guaranteed. (`#2883 <https://github.com/sunpy/sunpy/pull/2883>`__)
- Changed `sunpy.sun.models.interior` and `sunpy.sun.models.evolution` from `pandas.DataFrame` to `astropy.table.QTable` (`#2936 <https://github.com/sunpy/sunpy/pull/2936>`__)
- Minimum numpy version is now >=1.14.5 (`#2954 <https://github.com/sunpy/sunpy/pull/2954>`__)
- Removed ``sunpy.time.julian_day``, ``sunpy.time.julian_centuries``, ``sunpy.time.day_of_year``, ``sunpy.time.break_time``, ``sunpy.time.get_day``. (`#2999 <https://github.com/sunpy/sunpy/pull/2999>`__)
- Updated the solar values in `sunpy.sun.constants` to IAU 2015 values. (`#3001 <https://github.com/sunpy/sunpy/pull/3001>`__)
- Renamed `eccentricity_sunearth_orbit` to `eccentricity_sun_earth_orbit`. (`#3001 <https://github.com/sunpy/sunpy/pull/3001>`__)
- Renamed ``sunpy.image.rescale`` to `sunpy.image.resample`. (`#3044 <https://github.com/sunpy/sunpy/pull/3044>`__)
- Remove the ``basic_plot`` keyword argument from
  `~sunpy.map.Map.GenericMap.peek`. An example has been added to the gallery
  showing how to make a plot like this. (`#3109 <https://github.com/sunpy/sunpy/pull/3109>`__)
- `sunpy.map.GenericMap` will no longer use the key `solar_b0` as a value for heliographic latitude. (`#3115 <https://github.com/sunpy/sunpy/pull/3115>`__)
- `sunpy.map.GenericMap` now checks for a complete observer location rather than
  individually defaulting coordinates (lat, lon, distance) to Earth position. If
  any one of the three coordinates is missing from the header the observer will
  be defaulted to Earth and a warning raised. (`#3115 <https://github.com/sunpy/sunpy/pull/3115>`__)
- `sunpy.sun.sun` functions have been re-implemented using Astropy for significantly improved accuracy.  Some functions have been removed. (`#3137 <https://github.com/sunpy/sunpy/pull/3137>`__)
- All of the functions in `sunpy.sun.sun` and all of the Sun-specific functions in `sunpy.coordinates.ephemeris` have been moved to the new module `sunpy.coordinates.sun`. (`#3163 <https://github.com/sunpy/sunpy/pull/3163>`__)


Deprecations and Removals
-------------------------

- The deprecated ``sunpy.lightcurve``, ``sunpy.wcs`` and ``sunpy.spectra`` modules have now
  been removed. (`#2666 <https://github.com/sunpy/sunpy/pull/2666>`__)
- ``sunpy.instr.rhessi.get_obssumm_dbase_file`` ``sunpy.instr.rhessi.get_obssum_filename``, ``sunpy.instr.rhessi.get_obssumm_file`` have been removed. `Fido <sunpy.net.fido_factory.UnifiedDownloader>` should be used to download these files. (`#2808 <https://github.com/sunpy/sunpy/pull/2808>`__)
- Removed ``heliographic_solar_center`` in favour of `~sunpy.coordinates.ephemeris.get_sun_L0` and `~sunpy.coordinates.ephemeris.get_sun_B0` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``GenericClient.query`` in favour of `sunpy.net.dataretriever.GenericClient.search` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunearth_distance`` in favour of ``get_sunearth_distance`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``remove_lytaf_events_from_lightcurve`` in favour of `sunpy.instr.lyra.remove_lytaf_events_from_timeseries` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunpy.cm.get_cmap`` in favour of ``plt.get_cmap`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``database.query`` in favour of `sunpy.database.Database.search` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunpy.net.vso.InteractiveVSOClient`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``MapCube`` in favour of `~sunpy.map.MapSequence` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``solar_north`` in favour of ``get_sun_P`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``database.download`` in favour of `sunpy.database.Database.fetch` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunpy.map.GenericMap.pixel_to_data`` in favour of `sunpy.map.GenericMap.pixel_to_world` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``GenericClient.get`` in favour of `sunpy.net.dataretriever.GenericClient.fetch`. This changes applies to the other clients as well. (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed `Map.xrange` and `Map.yrange` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunpy.net.attrs.Wave`` in favour of `a.Wavelength <~sunpy.net.vso.attrs.Wavelength>` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``JSOCClient.check_request`` in favour of `drms.ExportRequest.status` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- `sunpy.net.vso.VSOClient.query_legacy` and `sunpy.net.vso.VSOClient.latest` have been deprecated as we strongly recommend people use `sunpy.net.Fido` for all queries. (`#2866 <https://github.com/sunpy/sunpy/pull/2866>`__)
- The deprecated ``sunpy.physics.transforms`` module has been removed, it is
  replaced by `sunpy.physics.solar_rotation` and
  `sunpy.physics.differential_rotation`. (`#2994 <https://github.com/sunpy/sunpy/pull/2994>`__)
- Removed `~sunpy.sun.sun.solar_cycle_number` because it was fundamentally flawed (`#3150 <https://github.com/sunpy/sunpy/pull/3150>`__)


Features
--------

- Change arguments to `sunpy.test` from ``offline=`` and ``online=`` to ``online`` and ``online_only``. This matches the behavior of the figure keyword arguments and comes as a part of a move to using a modified version of the Astropy test runner. (`#1983 <https://github.com/sunpy/sunpy/pull/1983>`__)
- asdf schemas and tags were added for the SunPy coordinate frames and `~sunpy.map.GenericMap` allowing these objects to be saved to and restored from `asdf <https://asdf.readthedocs.io/>`__ files. (`#2366 <https://github.com/sunpy/sunpy/pull/2366>`__)
- The images from image tests are now saved in a local folder for easy access. (`#2507 <https://github.com/sunpy/sunpy/pull/2507>`__)
- ``sunpy.map.MapCube`` has been renamed to `sunpy.map.MapSequence` to better reflect its use as a collection of map objects. (`#2603 <https://github.com/sunpy/sunpy/pull/2603>`__)
- Net search attributes now support tab completion of values and display a table of possible values when printed, to allow easier discoverability of possible search values. (`#2663 <https://github.com/sunpy/sunpy/pull/2663>`__)
- Running the figure tests now creates a page showing the differences between
  the expected figures and the figures produced from running the tests. (`#2681 <https://github.com/sunpy/sunpy/pull/2681>`__)
- Add support for Dask arrays in `sunpy.map.Map`. The map factory now checks a whitelist
  of array types rather than strictly checking if the array is of type `numpy.ndarray`. (`#2689 <https://github.com/sunpy/sunpy/pull/2689>`__)
- Persist the name of a coordinate, i.e. "earth" even though a concrete
  coordinate object has been calculated and use this string representation to change
  the way the sunpy frames are printed. This is primarily to facilitate displaying
  the name of the body rather than the concrete coordinate when printing a
  `~astropy.coordinates.SkyCoord`. (`#2723 <https://github.com/sunpy/sunpy/pull/2723>`__)
- `~sunpy.net.hek.HEKClient.search` now returns an `astropy.table.Table` instead of list of a `dict`. (`#2759 <https://github.com/sunpy/sunpy/pull/2759>`__)
- Add a downscaled HMI image to the sample data. (`#2782 <https://github.com/sunpy/sunpy/pull/2782>`__)
- Now able to create a `sunpy.map.Map` using an array and a `astropy.wcs.WCS` object. (`#2793 <https://github.com/sunpy/sunpy/pull/2793>`__)
- The download manager for `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloader.fetch>` has been replaced with
  `parfive <https://parfive.readthedocs.io/en/latest/>`__. This provides advanced
  progress bars, proper handling of overwriting and the ability to retry failed
  downloads. (`#2797 <https://github.com/sunpy/sunpy/pull/2797>`__)
- `sunpy.map.GenericMap` can now save out rice compressed FITS files. (`#2826 <https://github.com/sunpy/sunpy/pull/2826>`__)
- Now any SunPyDeprecationWarnings will cause an error when using pytest. (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Added full Tox support for SunPy tests, documentation build and figure tests. (`#2839 <https://github.com/sunpy/sunpy/pull/2839>`__)
- Transition the `sunpy.net.vso.VSOClient` from using suds to `zeep <https://python-zeep.readthedocs.io/en/master/>`__ as the SOAP
  library. This is a more actively maintained library, and should provide better
  support for the VSOs https endpoints. This change should have no effect on the
  public API of the `sunpy.net.vso.VSOClient`. (`#2866 <https://github.com/sunpy/sunpy/pull/2866>`__)
- Provided access to the Helioviewer header information using `~sunpy.net.helioviewer.HelioviewerClient.get_jp2_header` function. (`#2904 <https://github.com/sunpy/sunpy/pull/2904>`__)
- Add a new WSDL URL and port to support SunPy use of VSO instance at SDAC. (`#2912 <https://github.com/sunpy/sunpy/pull/2912>`__)
- Add support for COSMO K-Coronograph (KCOR) FITS data. (`#2916 <https://github.com/sunpy/sunpy/pull/2916>`__)
- Add logger messaging system based on `~astropy.logger.AstropyLogger`, cleaned up all warnings, removed all print statements. (`#2980 <https://github.com/sunpy/sunpy/pull/2980>`__)
- The function `sunpy.image.coalignment.get_correlation_shifts` now issues an error when the number of dimensions
  are not correct instead of a warning and returning None. (`#2980 <https://github.com/sunpy/sunpy/pull/2980>`__)
- The default location of the sunpy sample data has changed to be in the platform
  specific data directory as provided by `appdirs <https://github.com/ActiveState/appdirs>`__. (`#2993 <https://github.com/sunpy/sunpy/pull/2993>`__)
- Add timeseries support for EVE/ESP level 1 data in `sunpy.timeseries.sources.eve` (`#3032 <https://github.com/sunpy/sunpy/pull/3032>`__)
- The default style for Map plots have changed to reflect the changes in Astropy
  3.2. (`#3054 <https://github.com/sunpy/sunpy/pull/3054>`__)
- `sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst` can now account for light travel time when computing the (apparent) body position, as long as the observer location is provided. (`#3055 <https://github.com/sunpy/sunpy/pull/3055>`__)
- Added a helper function (`sunpy.map.make_fitswcs_header`) that allows users to create a meta header for custom created `sunpy.map.GenericMap`. (`#3083 <https://github.com/sunpy/sunpy/pull/3083>`__)
- Map plotting now accepts the optional keyword `clip_interval` for specifying a percentile interval for clipping.  For example, if the interval (5%, 99%) is specified, the bounds of the z axis are chosen such that the lowest 5% of pixels and the highest 1% of pixels are excluded. (`#3100 <https://github.com/sunpy/sunpy/pull/3100>`__)
- The new function `~sunpy.coordinates.get_horizons_coord` enables querying JPL HORIZONS for the locations of a wide range of solar-system bodies, including spacecraft. (`#3113 <https://github.com/sunpy/sunpy/pull/3113>`__)


Bug Fixes
---------

- Fix the bug that prevented VSO queries for HMI data from downloading file
  without specifying ``a.Physobs``. (`#2621 <https://github.com/sunpy/sunpy/pull/2621>`__)
- Fix `sunpy.map.mapcube.MapCube.plot`. The code had not been updated to support the changes to the wcsaxes helper functions. (`#2627 <https://github.com/sunpy/sunpy/pull/2627>`__)
- Replace all use of the deprecated ``sunpy.cm.get_cmap`` with `matplotlib.pyplot.get_cmap` to prevent deprecation warnings being raised. (`#2635 <https://github.com/sunpy/sunpy/pull/2635>`__)
- Fix generation of the coordinate transformation graph with Astropy 3.1.dev (`#2636 <https://github.com/sunpy/sunpy/pull/2636>`__)
- Prevent helioviewer from erroring when downloading file to a directory that
  does not exist. It will now create the directory when required. (`#2642 <https://github.com/sunpy/sunpy/pull/2642>`__)
- Fix transformations into/out of Heliographic Stonyhurst frame when 
  the coordinate representation is Cartesian. (`#2646 <https://github.com/sunpy/sunpy/pull/2646>`__)
- Running the figure tests with ``setup.py test`` now saves the figures and the hashes to the same directory as setup.py. (`#2658 <https://github.com/sunpy/sunpy/pull/2658>`__)
- `sunpy.instr.fermi.met_to_utc` now returns the correct utc time which takes into account the leap seconds that have passed. (`#2679 <https://github.com/sunpy/sunpy/pull/2679>`__)
- Support passing Python file objects to `sunpy.io.fits.write`. (`#2688 <https://github.com/sunpy/sunpy/pull/2688>`__)
- Added DRMS to setup.py so sunpy[all] installs it as a dependancy. (`#2693 <https://github.com/sunpy/sunpy/pull/2693>`__)
- Fix eve 0cs timeseries seperator regex to support Python 3.7 (`#2697 <https://github.com/sunpy/sunpy/pull/2697>`__)
- Fix the bug which crashes `~sunpy.map.sources.LASCOMap` for when 'date-obs' is reformatted agian from a self applied function. (`#2700 <https://github.com/sunpy/sunpy/pull/2700>`__)
- Change all instances of quantity_allclose to `astropy.units.allclose` this prevents pytest being needed to import `sunpy.coordinates` on Astropy 3 (`#2701 <https://github.com/sunpy/sunpy/pull/2701>`__)
- Fix RHESSI obssum file downloading to include the final day in the time range. (`#2714 <https://github.com/sunpy/sunpy/pull/2714>`__)
- Raise an error when transforming between HPC and HCC frames if the observer is not the same. (`#2725 <https://github.com/sunpy/sunpy/pull/2725>`__)
- Replaces the existing LASCO C2 and C3 color maps with new ones that perform better with JP2 and Level 0.5, 1 data. (`#2731 <https://github.com/sunpy/sunpy/pull/2731>`__)
- Do not attempt to save a FITS header comment for a keyword which is not in the header. This prevents an error on saving some maps after the metadata had been modified but not the comments. (`#2748 <https://github.com/sunpy/sunpy/pull/2748>`__)
- Add support for `~sunpy.map.sources.HMIMap` objects as input to `sunpy.instr.aia.aiaprep`. (`#2749 <https://github.com/sunpy/sunpy/pull/2749>`__)
- User can convert between HPC and HCC coordinates with different observers. This is implemented by automatically transforming the coordinate into HGS and then changing observer, and then transforming back to HCC. (`#2754 <https://github.com/sunpy/sunpy/pull/2754>`__)
- Changed default file type for Helioviewer to prevent decode errors. (`#2771 <https://github.com/sunpy/sunpy/pull/2771>`__)
- Increase figure size to avoid cutting off longer colormap names in `sunpy.cm.show_colormaps`. (`#2824 <https://github.com/sunpy/sunpy/pull/2824>`__)
- The sample data directory will no longer be created until files are downloaded
  to it. (`#2836 <https://github.com/sunpy/sunpy/pull/2836>`__)
- Timeseries and lightcurve will now respect updated config values for download directory. (`#2844 <https://github.com/sunpy/sunpy/pull/2844>`__)
- Always use _default_wrap_angle rather than hard coding a wrap angle in the init
  of a sunpy coordinate frame (`#2853 <https://github.com/sunpy/sunpy/pull/2853>`__)
- Ensure imageanimators only slice arrays with integers (`#2856 <https://github.com/sunpy/sunpy/pull/2856>`__)
- Fixed `sunpy.io.fits.write` to handle the keyword ``COMMENT`` correctly. (`#2880 <https://github.com/sunpy/sunpy/pull/2880>`__)
- If Carrington longitude ("crln_obs") is found in the FITS header, `~sunpy.map.Map` converts this to the correct Heliographic longitude. (`#2946 <https://github.com/sunpy/sunpy/pull/2946>`__)
- `sunpy.net.helio.hec.HECClient.time_query` now resolves the correct input time format. (`#2969 <https://github.com/sunpy/sunpy/pull/2969>`__)
- Fixes the calculation of the solar rotation of coordinates and the differential rotation of `sunpy.map.GenericMap`. (`#2972 <https://github.com/sunpy/sunpy/pull/2972>`__)
- Added back the FERMI GBM client to `sunpy.net.dataretriever.sources`. (`#2983 <https://github.com/sunpy/sunpy/pull/2983>`__)
- Fix bug in `sunpy.net.hek` which raised and error if a search returned zero results, now returns an empty `sunpy.net.hek.HEKTable`. (`#3046 <https://github.com/sunpy/sunpy/pull/3046>`__)
- `~sunpy.map.sources.AIAMap` now uses the provided HAE coordinates instead of the provided HGS coordinates to determine the observer location. (`#3056 <https://github.com/sunpy/sunpy/pull/3056>`__)
- Correctly zero pad milliseconds in the `sunpy.util.scraper.Scraper` formatting to prevent errors when the millisecond value was less than 100. (`#3063 <https://github.com/sunpy/sunpy/pull/3063>`__)
- Fix `sunpy.util.scraper.Scraper` failing if a directory is not found on a remote server. (`#3063 <https://github.com/sunpy/sunpy/pull/3063>`__)
- Correctly extract observer location from MDI and EIT data (`#3067 <https://github.com/sunpy/sunpy/pull/3067>`__)
- Fix HGS <> HCRS test due to Ecliptic frame changes in astropy 3.2 (`#3075 <https://github.com/sunpy/sunpy/pull/3075>`__)
- Fixes bug when creating a timeseries from a URL and bug when creating a TimeSeries from  older GOES/XRS fits files. (`#3081 <https://github.com/sunpy/sunpy/pull/3081>`__)
- Added `~sunpy.map.EUVIMap.rsun_obs`. It returns a quantity in arcsec consistent with other `sunpy.map.GenericMap` and overwrites mapbase's assumption of a photospheric limb as seen from Earth. (`#3099 <https://github.com/sunpy/sunpy/pull/3099>`__)
- Fixed bugs related to using `~sunpy.map.GenericMap.plot` and `~sunpy.map.GenericMap.peek` with the ``inline`` Matplotlib backend in Jupyter notebook. (`#3103 <https://github.com/sunpy/sunpy/pull/3103>`__)
- Make a correction to `sunpy.coordinates.wcs_utils.solar_wcs_frame_mapping` so
  that `astropy.wcs.WCS` objects are correctly converted to
  `sunpy.coordinates.frames` objects irrespective of the ordering of the axes. (`#3116 <https://github.com/sunpy/sunpy/pull/3116>`__)
- The `solar_rotate_coordinate` function returns a coordinate that accounts for the location of the new observer. (`#3123 <https://github.com/sunpy/sunpy/pull/3123>`__)
- Add support for rotation parameters to `sunpy.map.make_fitswcs_header`. (`#3139 <https://github.com/sunpy/sunpy/pull/3139>`__)
- Improve the implementation of `~sunpy.physics.differential_rotation.differential_rotate` the image warping when transforming Maps for differential rotation and change in observer position. (`#3149 <https://github.com/sunpy/sunpy/pull/3149>`__)
- Fix a bug where new helioviewer sources potentially cause `~sunpy.net.helioviewer.HelioviewerClient.data_sources` to error. (`#3162 <https://github.com/sunpy/sunpy/pull/3162>`__)


Improved Documentation
----------------------

- Organise the gallery into sections based on example type and tidy up a little. (`#2624 <https://github.com/sunpy/sunpy/pull/2624>`__)
- Added gallery example showing the conversion of Helioprojective Coordinates to Altitude/Azimuth Coordinates to and back. (`#2656 <https://github.com/sunpy/sunpy/pull/2656>`__)
- Add contribution guidelines for the sunpy example gallery. (`#2682 <https://github.com/sunpy/sunpy/pull/2682>`__)
- Added a gallery example for "Downloading and plotting a HMI image" and "Creating a Composite map". (`#2746 <https://github.com/sunpy/sunpy/pull/2746>`__)
- Added an example for `~sunpy.visualization.animator.ImageAnimatorWCS`. (`#2752 <https://github.com/sunpy/sunpy/pull/2752>`__)
- Minor changes to the developer guide regarding sprint labels. (`#2765 <https://github.com/sunpy/sunpy/pull/2765>`__)
- Copyedited and corrected the solar cycles example. (`#2770 <https://github.com/sunpy/sunpy/pull/2770>`__)
- Changed "online" mark to "remote_data" and made formatting of marks consistent. (`#2799 <https://github.com/sunpy/sunpy/pull/2799>`__)
- Add a missing plot to the end of the units and coordinates guide. (`#2813 <https://github.com/sunpy/sunpy/pull/2813>`__)
- Added gallery example showing how to access the SunPy colormaps (`#2865 <https://github.com/sunpy/sunpy/pull/2865>`__)
- Added gallery example showing how to access the SunPy solar physics constants. (`#2882 <https://github.com/sunpy/sunpy/pull/2882>`__)
- Major clean up of the developer documentation. (`#2951 <https://github.com/sunpy/sunpy/pull/2951>`__)
- Overhaul of the install intructions for the guide section of our documentation. (`#3147 <https://github.com/sunpy/sunpy/pull/3147>`__)


Trivial/Internal Changes
------------------------

- `~sunpy.time.parse_time` now uses `singledispatch` underneath. (`#2408 <https://github.com/sunpy/sunpy/pull/2408>`__)
- Revert the handling of ``quantity_allclose`` now that `astropy/astropy#7252 <https://github.com/astropy/astropy/pull/7252>`__ is merged. This also bumps the minimum astropy version to 3.0.2. (`#2598 <https://github.com/sunpy/sunpy/pull/2598>`__)
- Replace the subclasses of matplotlib Slider and Button in `sunpy.visualization` with partial functions. (`#2613 <https://github.com/sunpy/sunpy/pull/2613>`__)
- Sort the ana C source files before building to enable reproducible builds. (`#2637 <https://github.com/sunpy/sunpy/pull/2637>`__)
- We are now using `towncrier <https://github.com/hawkowl/towncrier>`__ to
  generate our changelogs. (`#2644 <https://github.com/sunpy/sunpy/pull/2644>`__)
- Moved figure tests to Python 3.6. (`#2655 <https://github.com/sunpy/sunpy/pull/2655>`__)
- Removed old metaclass used for Map and TimeSeries as we have now moved to Python 3.6. (`#2655 <https://github.com/sunpy/sunpy/pull/2655>`__)
- Updated astropy_helpers to v3.0.2. (`#2655 <https://github.com/sunpy/sunpy/pull/2655>`__)
- When running image tests, a comparison HTML page is now generated to show
  the generated images and expected images. (`#2660 <https://github.com/sunpy/sunpy/pull/2660>`__)
- Change to using pytest-cov for coverage report generation to enable support for parallel builds (`#2667 <https://github.com/sunpy/sunpy/pull/2667>`__)
- Use of `textwrap` to keep source code indented when multiline texts is used (`#2671 <https://github.com/sunpy/sunpy/pull/2671>`__)
- Fix mispelling of private attribute ``_default_heliographic_latitude`` in map. (`#2730 <https://github.com/sunpy/sunpy/pull/2730>`__)
- Miscellaneous fixes to developer docs about building sunpy's documentation. (`#2825 <https://github.com/sunpy/sunpy/pull/2825>`__)
- Changed `sunpy.instr.aia.aiaprep` to update BITPIX keyword to reflect the float64 dtype. (`#2831 <https://github.com/sunpy/sunpy/pull/2831>`__)
- Remove warning from ``GenericMap.submap`` when using pixel ``Quantities`` as input. (`#2833 <https://github.com/sunpy/sunpy/pull/2833>`__)
- Remove the usage of six and all ``__future__`` imports (`#2837 <https://github.com/sunpy/sunpy/pull/2837>`__)
- Fix SunPy Coordinate tests with Astropy 3.1 (`#2838 <https://github.com/sunpy/sunpy/pull/2838>`__)
- Stores entries from directories into database sorted by name. It adds mocks to the database user guide examples. (`#2873 <https://github.com/sunpy/sunpy/pull/2873>`__)
- Fix all DeprecationWarning: invalid escape sequence. (`#2885 <https://github.com/sunpy/sunpy/pull/2885>`__)
- Used `unittest.mock` for creating offline tests for simulating online tests for `test_noaa.py` (`#2900 <https://github.com/sunpy/sunpy/pull/2900>`__)
- Fix support for pip 19 and isolated builds (`#2915 <https://github.com/sunpy/sunpy/pull/2915>`__)
- Moved to using `AppDirs <https://github.com/ActiveState/appdirs>`__ as the place to host our configuration file. (`#2922 <https://github.com/sunpy/sunpy/pull/2922>`__)
- Users can now use fewer keywords in our `~sunpy.net.HelioviewerClient` to access the available sources. Either by `observatory` and `measurement` or `instrument` and `measurement` as this much information is enough to get the source ID for most of the cases. (`#2926 <https://github.com/sunpy/sunpy/pull/2926>`__)
- Remove the pytest dependancy on the ``GenericMap`` asdf tag. (`#2943 <https://github.com/sunpy/sunpy/pull/2943>`__)
- Fix initialization of `~sunpy.net.vso.VSOClient` when no WSDL link is found. (`#2981 <https://github.com/sunpy/sunpy/pull/2981>`__)


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

- `sunpy.coordinates.representation` has been removed. Longitude wrapping is now done in the constructor of the frames. [#2431]
- Propagation of ``obstime`` in the coordinate frame transformation has changed, this means in general when transforming directly between frames (not `~astropy.coordinates.SkyCoord`) you will have to specify ``obstime`` in more places. [#2461]
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
   https://docs.astropy.org/en/latest/constants/index.html
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
