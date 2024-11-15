6.0.0 (2024-07-19)
==================

Breaking Changes
----------------

- Arguments for :meth:`~sunpy.map.GenericMap.reproject_to` after the target WCS are now keyword-only. (`#7339 <https://github.com/sunpy/sunpy/pull/7339>`__)
- Arguments for :meth:`sunpy.timeseries.GenericTimeSeries.peek` are now keywords only. (`#7340 <https://github.com/sunpy/sunpy/pull/7340>`__)
- Removed scikit-image from the "image" extra group and created a new "scikit-image" extra group. (`#7536 <https://github.com/sunpy/sunpy/pull/7536>`__)
- The "all" extra group now will install all optional packages.

  This now includes the following packages:

  - asdf
  - glmyur
  - opencv
  - scikit-image
  - spiceypy (`#7536 <https://github.com/sunpy/sunpy/pull/7536>`__)
- Removed the "dask" extra group. (`#7536 <https://github.com/sunpy/sunpy/pull/7536>`__)
- ``sunpy.io.read_file`` and ``sunpy.io.write_file`` are deprecated and will be removed in the future.
  These were intended to be private functions and should not be used. (`#7537 <https://github.com/sunpy/sunpy/pull/7537>`__)
- The ANA C code has been deprecated (`sunpy.io.ana.read`, `sunpy.io.ana.get_header`, `sunpy.io.ana.write`) and may be removed in a future sunpy release.
  Please contact us here: https://community.openastronomy.org/t/possible-deprecation-of-ana-file-readers-and-writers-in-sunpy if you are making use of this code. (`#7642 <https://github.com/sunpy/sunpy/pull/7642>`__)
- The `.EUIMap` class now returns the ``DATE-BEG`` key for `.GenericMap.date` while continuing to use ``DATE-AVG`` as the reference date for the coordinate system. (`#7682 <https://github.com/sunpy/sunpy/pull/7682>`__)
- The `.GenericMap.date` key priority order has changed to be consistent with it representing the "canonical" observation time.
  ``DATE-OBS`` continues to have the highest priority, but now ``DATE-BEG`` has higher priority than ``DATE-AVG``. (`#7682 <https://github.com/sunpy/sunpy/pull/7682>`__)
- A new property `.GenericMap.reference_date` has been added to decouple the reference date for the coordinate system from the "canonical" observation time.
  This new property is now passed through to the map's WCS object as ``dateavg`` and is the time used for `.GenericMap.coordinate_frame` and `.GenericMap.observer_coordinate`. (`#7682 <https://github.com/sunpy/sunpy/pull/7682>`__)


Deprecations
------------

- :meth:`~sunpy.coordinates.Helioprojective.assume_spherical_screen` has been deprecated in favor of `~sunpy.coordinates.SphericalScreen`. (`#7115 <https://github.com/sunpy/sunpy/pull/7115>`__)
- :func:`sunpy.physics.differential_rotation.diff_rot` has been deprecated and replaced by :func:`sunpy.sun.models.differential_rotation`. (`#7409 <https://github.com/sunpy/sunpy/pull/7409>`__)
- Deprecated all positional arguments in :meth:`sunpy.map.GenericMap.plot` method.
  The ``annotate``, ``axes``, ``title``, ``clip_interval`` arguments should be passed as keyword arguments (e.g., ``..., title=True, ...``) instead. (`#7421 <https://github.com/sunpy/sunpy/pull/7421>`__)
- The keyword ``response_format`` in :meth:`sunpy.net.vso.VSOClient.search` has been deprecated.
  This was introduced to preserve legacy behaviour of the VSO client, to return
  ``sunpy.net.vso.legacy_response.QueryResponse`` instead of `sunpy.net.vso.table_response.VSOQueryResponseTable` objects.
  This behaviour has been the default for over 4 years and the keyword is no longer needed.
  This keyword and the older ``sunpy.net.vso.legacy_response.QueryResponse`` class will be removed in sunpy 7.0.
  The keyword ``progress`` in :meth:`sunpy.net.hek2vso.H2VClient.full_query` has been deprecated and will be removed in sunpy 7.0. (`#7468 <https://github.com/sunpy/sunpy/pull/7468>`__)


Removals
--------

- ``sunpy.database`` has been removed. (`#7320 <https://github.com/sunpy/sunpy/pull/7320>`__)
- ``sunpy.map.header_helper.meta_keywords`` has been removed. (`#7337 <https://github.com/sunpy/sunpy/pull/7337>`__)
- ``sunpy.net.helioviewer.HelioviewerClient`` has been removed. Use the `hvpy <https://hvpy.readthedocs.io/en/latest/>`__ package instead. (`#7338 <https://github.com/sunpy/sunpy/pull/7338>`__)
- There was a private "Maxwell" unit within `sunpy.map` to register it before astropy had support for it.
  This has now been removed in favour of using the astropy version. (`#7383 <https://github.com/sunpy/sunpy/pull/7383>`__)


New Features
------------

- ``sunpy.io.read_file`` will now try to detect the filetype based on the content and then fallback to using the file extension. (`#6736 <https://github.com/sunpy/sunpy/pull/6736>`__)
- It is now possible to read the comments in a header from a JPEG2000 file. (`#6841 <https://github.com/sunpy/sunpy/pull/6841>`__)
- Added the ability for `sunpy.map.Map` to load files from a generator. (`#7024 <https://github.com/sunpy/sunpy/pull/7024>`__)
- Added `~sunpy.coordinates.PlanarScreen` for interpreting 2D `~sunpy.coordinates.Helioprojective` coordinates as being on the inside of a planar screen. (`#7115 <https://github.com/sunpy/sunpy/pull/7115>`__)
- Added the ability to pass ``clip_interval`` to :meth:`sunpy.map.mapsequence.MapSequence.plot`. (`#7253 <https://github.com/sunpy/sunpy/pull/7253>`__)
- Add support for the ``fill`` keyword in :meth:`~sunpy.map.GenericMap.draw_contours` to allow for filled contours. (`#7281 <https://github.com/sunpy/sunpy/pull/7281>`__)
- :func:`~sunpy.coordinates.get_horizons_coord` now supports time arrays with up to 10,000 elements. (`#7319 <https://github.com/sunpy/sunpy/pull/7319>`__)
- Add an example of plotting a rectangle on a map with a rotation angle relative to the axes (:ref:`sphx_glr_generated_gallery_plotting_plot_rotated_rectangle.py`). (`#7348 <https://github.com/sunpy/sunpy/pull/7348>`__)
- Added testing and explicit support for Python 3.12. (`#7351 <https://github.com/sunpy/sunpy/pull/7351>`__)
- Added warning when importing a submodule without installing that submodules extra dependencies. (`#7369 <https://github.com/sunpy/sunpy/pull/7369>`__)
- Added a warning message for ``rsun`` mismatch in :meth:`~sunpy.map.GenericMap.reproject_to` method. (`#7370 <https://github.com/sunpy/sunpy/pull/7370>`__)
- Added a new optional extra group to install "opencv" if you want to it for affine transforms.

  .. code-block:: bash

      pip install sunpy[opencv] (`#7383 <https://github.com/sunpy/sunpy/pull/7383>`__)
- Increased minimum versions for:

  - asdf >= 2.12.0
  - asdf-astropy >= 0.2.0
  - astropy >= 5.2.0
  - beautifulsoup4 >= 4.11.0
  - cdflib >= 0.4.4
  - dask >= 2022.5.2
  - h5netcdf > =1.0.0
  - h5py >= 3.7.0
  - lxml >= 4.9.0
  - opencv-python >= 4.6.0.66
  - pandas >= 1.4.0
  - python >= 3.10
  - reproject >= 0.9.0
  - requests >= 2.28.0
  - scikit-image >= 0.19.0
  - scipy >= 1.8.0
  - spiceypy >= 5.0.0
  - tqdm >= 4.64.0
  - zeep >= 4.1.0 (`#7383 <https://github.com/sunpy/sunpy/pull/7383>`__)
- :meth:`sunpy.map.GenericMap.draw_contours` don't run internal transform code if ``transform`` keyword is provided. (`#7427 <https://github.com/sunpy/sunpy/pull/7427>`__)
- Update ASDF schemas for upcoming ASDF standard 1.6.0. (`#7432 <https://github.com/sunpy/sunpy/pull/7432>`__)
- Add a new map source `~sunpy.map.sources.gong.GONGHalphaMap` for GONG H-Alpha data. (`#7451 <https://github.com/sunpy/sunpy/pull/7451>`__)
- Added :func:`~sunpy.coordinates.spice.get_rotation_matrix` to obtain the rotation matrix between the orientations of two SPICE frames, which is particularly useful for transforming vector fields. (`#7452 <https://github.com/sunpy/sunpy/pull/7452>`__)
- Allow units to be passed to `~sunpy.map.header_helper.make_fitswcs_header` as strings. (`#7454 <https://github.com/sunpy/sunpy/pull/7454>`__)
- A new client (`sunpy.net.dataretriever.ADAPTClient`) has been added to search and download `ADAPT <https://gong.nso.edu/adapt/maps/gong/>`__ files. (`#7463 <https://github.com/sunpy/sunpy/pull/7463>`__)
- `sunpy.net.jsoc.JSOCClient` queries now return the SUMS directory paths as the segment key value in the results table. (`#7469 <https://github.com/sunpy/sunpy/pull/7469>`__)
- Allow the screen radius to be set when using `~sunpy.coordinates.SphericalScreen`. (`#7532 <https://github.com/sunpy/sunpy/pull/7532>`__)
- Added a "core" extra group that does not install any truly optional dependencies.
  It only includes the dependencies that are required to import sunpy and all subpackages.

  This means it will not install:

  - asdf
  - glymur
  - opencv
  - scikit-image
  - spiceypy (`#7536 <https://github.com/sunpy/sunpy/pull/7536>`__)
- Updated :meth:`sunpy.map.GenericMap.submap` to check if it is about to work on locations with NaNs now errors and informs the user that they likely want to use :meth:`~sunpy.coordinates.Helioprojective.assume_spherical_screen` so that the off-disk 2D coordinate can be converted to a 3D coordinate. (`#7543 <https://github.com/sunpy/sunpy/pull/7543>`__)
- `~sunpy.map.GenericMap` will now assign units of DN without a warning or error. (`#7585 <https://github.com/sunpy/sunpy/pull/7585>`__)
- Add a new map source `~sunpy.map.sources.ADAPTMap` for ADvanced Adaptive Prediction Technique (ADAPT) data files. (`#7640 <https://github.com/sunpy/sunpy/pull/7640>`__)
- Added support for JSOC's HMI millisecond TAI time format.
  Previously, it would only work with seconds. (`#7656 <https://github.com/sunpy/sunpy/pull/7656>`__)
- Added build support for aarch64 wheels. (`#7679 <https://github.com/sunpy/sunpy/pull/7679>`__)


Bug Fixes
---------

- Long object names are no longer truncated in the logging output of :func:`~sunpy.coordinates.get_horizons_coord`. (`#7319 <https://github.com/sunpy/sunpy/pull/7319>`__)
- When calling :meth:`sunpy.map.GenericMap.rotate` on an integer data array, with ``missing`` set to NaN (the default value), the method will now itself raise an informative error message instead deferring to NumPy to raise the error. (`#7344 <https://github.com/sunpy/sunpy/pull/7344>`__)
- Fixed the appearance of a double "Notes" heading in `~sunpy.map.Map` subclasses. (`#7376 <https://github.com/sunpy/sunpy/pull/7376>`__)
- `~sunpy.map.Map` with UINT8 data will now not error on plotting due to normalization.
  We now skip adding a normalization. (`#7422 <https://github.com/sunpy/sunpy/pull/7422>`__)
- When calling :meth:`~sunpy.map.GenericMap.reproject_to` along with both context managers :func:`~sunpy.coordinates.propagate_with_solar_surface` and :meth:`~sunpy.coordinates.Helioprojective.assume_spherical_screen` now raises a warning. (`#7437 <https://github.com/sunpy/sunpy/pull/7437>`__)
- Fix a bug which caused ``Fido.search`` to crash due to SSL certificate verification error for the `~sunpy.net.helio.HECClient` now returns no results and logs a warning in this case. (`#7446 <https://github.com/sunpy/sunpy/pull/7446>`__)
- Fixed the sanitization of the names of files downloaded via VSO so that periods are no longer replaced and case is no longer forced to be lowercase. (`#7453 <https://github.com/sunpy/sunpy/pull/7453>`__)
- The creation of the series string for a JSOC query was not adding the correct escape characters for  comparison values for keywords.
  This was causing the JSOC to error. (`#7467 <https://github.com/sunpy/sunpy/pull/7467>`__)
- The EVE L0CS client now uses the new URLs for the data from LASP. (`#7483 <https://github.com/sunpy/sunpy/pull/7483>`__)
- JPEG2000 files are now saved with the correct orientation. Previously they would be vertically flipped when saved. (`#7486 <https://github.com/sunpy/sunpy/pull/7486>`__)
- Fixed a very minor inaccuracy in three `sunpy.map` utility functions (:func:`~sunpy.map.contains_full_disk`, :func:`~sunpy.map.coordinate_is_on_solar_disk`, and :func:`~sunpy.map.is_all_off_disk`) resulting from the accidental use of the small-angle approximation. (`#7512 <https://github.com/sunpy/sunpy/pull/7512>`__)
- The :meth:`~sunpy.map.GenericMap.rotate` function now correctly updates the NAXISi. (`#7522 <https://github.com/sunpy/sunpy/pull/7522>`__)
- Added a check in `sunpy.physics.differential_rotation.solar_rotate_coordinate` to ensure the input frame has an "observer" attribute before replicating frame
  attributes, preventing potential issues with frames lacking this attribute. (`#7526 <https://github.com/sunpy/sunpy/pull/7526>`__)
- Fixed an inaccuracy in the implementation of `~sunpy.coordinates.HeliocentricEarthEcliptic` and `~sunpy.coordinates.GeocentricSolarEcliptic` such that the Earth was not exactly in the XY plane, but rather had an error of up ~10 meters. (`#7530 <https://github.com/sunpy/sunpy/pull/7530>`__)
- The maximum records in `~sunpy.net.helio.HECClient` now are 20000. (`#7540 <https://github.com/sunpy/sunpy/pull/7540>`__)
- Fixed a bug with any coordinate transformation starting in `~sunpy.coordinates.GeocentricEarthEquatorial` (GEI) returning output with AU as the length unit, rather than preserving the length unit of the initial coordinate. (`#7545 <https://github.com/sunpy/sunpy/pull/7545>`__)
- Fixed a bug that interfered with :func:`astropy.wcs.utils.celestial_frame_to_wcs` when working with a custom subclass of :class:`~sunpy.coordinates.frames.SunPyBaseCoordinateFrame`. (`#7594 <https://github.com/sunpy/sunpy/pull/7594>`__)
- Fixed bug where conversion of results from the HEKClient to Astropy Time failed when some values where empty or missing for the values of event_strattime, event_endtime or event_peaktime (`#7627 <https://github.com/sunpy/sunpy/pull/7627>`__)
- Fix the `~sunpy.map.sources.gong.GONGHalphaMap.rsun_obs` to use correct header information ``solar-r`` keyword. (`#7652 <https://github.com/sunpy/sunpy/pull/7652>`__)
- Fix compilation with gcc 14, avoid implicit pointer conversions. (`#7662 <https://github.com/sunpy/sunpy/pull/7662>`__)
- Fixed a bug where "DN" was not able to be parsed by `~sunpy.map.header_helper.make_fitswcs_header` due to strict checking
  against the FITS standard. This is now consistent with how unit strings are parsed in `~sunpy.map.GenericMap`. (`#7730 <https://github.com/sunpy/sunpy/pull/7730>`__)
- Fixed a bug where `~sunpy.map.sources.XRTMap` was still defaulting to counts rather than DN. (`#7744 <https://github.com/sunpy/sunpy/pull/7744>`__)


Documentation
-------------

- Added a how-to guide for manipulating grid lines on `~sunpy.map.GenericMap`. (`#6978 <https://github.com/sunpy/sunpy/pull/6978>`__)
- Created a how to guide on fixing metadata that is either missing or incorrect before passing the header into the `~sunpy.map.Map` class. (`#7262 <https://github.com/sunpy/sunpy/pull/7262>`__)
- Fixed the usage of :meth:`~sunpy.map.GenericMap.superpixel` in :ref:`sphx_glr_generated_gallery_map_map_resampling_and_superpixels.py`. (`#7316 <https://github.com/sunpy/sunpy/pull/7316>`__)
- Added Clarification on setting JSOC Email. (`#7329 <https://github.com/sunpy/sunpy/pull/7329>`__)
- Added explanation text to :ref:`sphx_glr_generated_gallery_plotting_plotting_blank_map.py` about the offset between "(0, 0)" in helioprojective coordinates and the heliographic equator. (`#7352 <https://github.com/sunpy/sunpy/pull/7352>`__)
- Convert draw rectangle gallery example into a how-to guide(:ref:`sunpy-how-to-create-rectangle-on-map`) (`#7435 <https://github.com/sunpy/sunpy/pull/7435>`__)
- Fix a VSO doctest due to VSO now returning level one EIT data. (`#7483 <https://github.com/sunpy/sunpy/pull/7483>`__)
- Add an example gallery entry demonstrating how to use the coordinates framework to compute intersections
  between instrument lines of sight and a simulation domain. (`#7491 <https://github.com/sunpy/sunpy/pull/7491>`__)
- Updated the examples for :func:`~sunpy.visualization.colormaps.color_tables.hmi_mag_color_table` that used older styles of plotting (`#7692 <https://github.com/sunpy/sunpy/pull/7692>`__)


Internal Changes
----------------

- :meth:`sunpy.net.jsoc.JSOCClient.fetch` called ``drms`` API that passed a ``progress`` keyword which added extra print statements to the console.
  This has been removed in ``drms`` 0.7.0, which had breaking API changes within this release.
  As a result, we increased the minimum required version of ``drms`` to 0.7.1.

  This specifically refers to the following information that was printed to the console by default:

  ``"Export request pending. [id=X, status=X]"``
  ``"Waiting for X seconds..."``
  ``"Request not found on server, X retries left."``

  These were handled by ``drms`` and are now logging messages.

  If you want to silence these messages, you can set the logging level to ``WARNING`` or higher.

  .. code-block:: python

      import logging
      drms_logger = logging.getLogger("drms")
      drms_logger.setLevel(logging.WARNING)

      from sunpy.net import fido, attrs

  Note, you have to do it before you import ``fido``. (`#7307 <https://github.com/sunpy/sunpy/pull/7307>`__)
- The function :func:`~sunpy.coordinates.get_horizons_coord` no longer calls the ``astroquery`` package, so ``astroquery`` is no longer a dependency. (`#7319 <https://github.com/sunpy/sunpy/pull/7319>`__)
- The ``requests`` package is a now formally a core dependency.
  ``requests`` was already commonly installed as an implied dependency of `sunpy.net` or for building documentation. (`#7319 <https://github.com/sunpy/sunpy/pull/7319>`__)
- `~sunpy.net.jsoc.attrs.Notify` checks that a valid email address has been given as a value. (`#7342 <https://github.com/sunpy/sunpy/pull/7342>`__)
- The ``delim_whitespace`` keyword in `pandas.read_csv` is deprecated and was updated with ``sep='\s+'``.
  This should have no affect on the output of the code. (`#7350 <https://github.com/sunpy/sunpy/pull/7350>`__)
- Fixed an environment-specific failure of a unit test for :meth:`sunpy.coordinates.Helioprojective.is_visible`. (`#7356 <https://github.com/sunpy/sunpy/pull/7356>`__)
- Moved to ``pyproject.toml`` and removed ``setup.py`` and ``setup.cfg``. (`#7384 <https://github.com/sunpy/sunpy/pull/7384>`__)
- ``pyerfa`` is now a new direct dependency.
  It has been an indirect dependency from sunpy 3.1, over two years ago. (`#7397 <https://github.com/sunpy/sunpy/pull/7397>`__)
- Increased Python minimum version to be >= 3.10. (`#7402 <https://github.com/sunpy/sunpy/pull/7402>`__)
- Fixed an unnecessary division computation when performing a unsupported division operation using a `~sunpy.map.Map`. (`#7551 <https://github.com/sunpy/sunpy/pull/7551>`__)
- Updated the internal URL for the `~sunpy.net.dataretriever.sources.norh.NoRHClient` to point to a HTTPS archive of the NoRH data. (`#7696 <https://github.com/sunpy/sunpy/pull/7696>`__)


5.1.0 (2023-11-20)
==================

New Features
------------

- Added the ability to skip over errors raised for invalid fits files when passing a list of files to map using the existing keyword argument ``silence_errors``. (`#7018 <https://github.com/sunpy/sunpy/pull/7018>`__)
- Added a :meth:`sunpy.coordinates.Helioprojective.is_visible` method to return whether the coordinate is visible (i.e., not obscured from the observer assuming that the Sun is an opaque sphere). (`#7118 <https://github.com/sunpy/sunpy/pull/7118>`__)
- Added a keyword option (``quiet``) for :func:`~sunpy.coordinates.get_body_heliographic_stonyhurst` to silence the normal reporting of the light-travel-time correction when ``observer`` is specified. (`#7142 <https://github.com/sunpy/sunpy/pull/7142>`__)
- Added the function :func:`sunpy.coordinates.sun.eclipse_amount` to calculate the solar-eclipse amount for an observer. (`#7142 <https://github.com/sunpy/sunpy/pull/7142>`__)
- Add a keyword (``map_center_longitude``) to :func:`~sunpy.map.header_helper.make_heliographic_header` for centering the heliographic map at a longitude other than zero longitude. (`#7143 <https://github.com/sunpy/sunpy/pull/7143>`__)
- The minimum required version of ``Glymur`` (an optional dependency for reading JPEG2000 files) has been increase to 0.9.1. (`#7164 <https://github.com/sunpy/sunpy/pull/7164>`__)
- Added new default colormap scalings for WISPR Maps. Plots are now clipped at zero, and `~astropy.visualization.AsinhStretch` is used for the scaling to ensure coronal details are visible despite the much-brighter stars. Parsing of the ``detector`` and ``level`` fields of the FITS headers is also improved. (`#7180 <https://github.com/sunpy/sunpy/pull/7180>`__)
- When creating a coordinate or coordinate frame without specifying ``obstime``, the ``obstime`` value from the ``observer`` frame attribute will be used if present. (`#7186 <https://github.com/sunpy/sunpy/pull/7186>`__)
- Added a GONG synoptic map class which fixes non-compliant FITS metadata (`#7220 <https://github.com/sunpy/sunpy/pull/7220>`__)
- Added the module `sunpy.coordinates.spice` to enable the use of the `~astropy.coordinates.SkyCoord` API to perform computations using `SPICE <https://naif.jpl.nasa.gov/naif/>`__ kernels. (`#7237 <https://github.com/sunpy/sunpy/pull/7237>`__)
- Added three coordinate frames that depend on the orientation of Earth's magnetic dipole: `~sunpy.coordinates.Geomagnetic` (MAG), `~sunpy.coordinates.SolarMagnetic` (SM), and `~sunpy.coordinates.GeocentricSolarMagnetospheric` (GSM). (`#7239 <https://github.com/sunpy/sunpy/pull/7239>`__)


Bug Fixes
---------

- Fix RHESSI (`~sunpy.net.dataretriever.RHESSIClient`) fallback server detection. (`#7092 <https://github.com/sunpy/sunpy/pull/7092>`__)
- Fix bug in :func:`~sunpy.coordinates.get_horizons_coord` when specifying a time range via a dictionary that could cause the returned times to be slightly different from the supplied times. (`#7106 <https://github.com/sunpy/sunpy/pull/7106>`__)
- Updated the url of the `~sunpy.net.dataretriever.GBMClient` to match on files other than those that end with version 0 (i.e., V0.pha). (`#7148 <https://github.com/sunpy/sunpy/pull/7148>`__)
- When directly instantiating a `~astropy.wcs.WCS` from a FITS header that contains both Stonyhurst and Carrington heliographic coordinates for the observer location, the Stonyhurst coordinates will now be prioritized.
  This behavior is now consistent with the `~sunpy.map.Map` class, which has always prioritized Stonyhurst coordinates over Carrington coordinates. (`#7188 <https://github.com/sunpy/sunpy/pull/7188>`__)
- Fixed a bug with :func:`~sunpy.map.sample_at_coords()` where sampling outside the bounds of the map would sometimes not error and instead return strange pixel values. (`#7206 <https://github.com/sunpy/sunpy/pull/7206>`__)
- Improved code when loading CDF files to improve performance and avoid raising of pandas performance warnings. (`#7247 <https://github.com/sunpy/sunpy/pull/7247>`__)
- Fixed a bug with :meth:`sunpy.map.GenericMap.plot` where setting ``norm`` to ``None`` would result in an error. (`#7261 <https://github.com/sunpy/sunpy/pull/7261>`__)


Documentation
-------------

- Removed the specification of a non-identity rotation matrix in two reprojection examples. (`#7114 <https://github.com/sunpy/sunpy/pull/7114>`__)
- Added an example (:ref:`sphx_glr_generated_gallery_showcase_stereoscopic_3d.py`) for how to make an anaglyph 3D (i.e., red-cyan) image from a stereoscopic observation. (`#7123 <https://github.com/sunpy/sunpy/pull/7123>`__)
- Added an example (:ref:`sphx_glr_generated_gallery_showcase_eclipse_amount.py`) to show how to obtain information about a solar eclipse using :func:`sunpy.coordinates.sun.eclipse_amount`. (`#7142 <https://github.com/sunpy/sunpy/pull/7142>`__)
- Changed the :ref:`sphx_glr_generated_gallery_map_masking_hmi.py` to reproject AIA to HMI instead of the other way around.
  This is to avoid interpolating the HMI LOS magnetic field data. (`#7160 <https://github.com/sunpy/sunpy/pull/7160>`__)
- Fixed the timeseries peak finding example.
  Previously there was a bug when plotting the data with pandas. (`#7199 <https://github.com/sunpy/sunpy/pull/7199>`__)
- Added an example (:ref:`sphx_glr_generated_gallery_units_and_coordinates_spice.py`) for how to perform `SPICE <https://naif.jpl.nasa.gov/naif/>`__ computations using the `~astropy.coordinates.SkyCoord` API. (`#7237 <https://github.com/sunpy/sunpy/pull/7237>`__)


Deprecations
------------

- Deprecated ``silence_errors`` in Map and Timeseries.
  This has been replaced with ``allow_errors`` keyword. (`#7021 <https://github.com/sunpy/sunpy/pull/7021>`__)
- The ``sunpy.coordinates.transformations`` module is now slated for removal from the public API as it consists of internal functions used by coordinate transformations.
  The context managers :func:`sunpy.coordinates.transform_with_sun_center` and :func:`sunpy.coordinates.propagate_with_solar_surface` should be accessed under `sunpy.coordinates`. (`#7113 <https://github.com/sunpy/sunpy/pull/7113>`__)


Removals
--------

- ``sunpy.map.extract_along_coord()`` has been removed.
  Instead, use :func:`~sunpy.map.pixelate_coord_path`, and then pass its output to :func:`~sunpy.map.sample_at_coords`.
  ``pixelate_coord_path`` uses a different line algorithm by default, but you can specify ``bresenham=True`` as an argument to use the same line algorithm as ``extract_along_coord``. (`#7200 <https://github.com/sunpy/sunpy/pull/7200>`__)
- ``sunpy.visualisation.limb.draw_limb()`` has been removed.
  Use :func:`sunpy.visualization.drawing.limb` instead. (`#7202 <https://github.com/sunpy/sunpy/pull/7202>`__)
- Removed ``GenericTimeSeries.index``.
  Use ``GenericTimeseries.time`` instead as a direct replacement. (`#7203 <https://github.com/sunpy/sunpy/pull/7203>`__)
- Removed the deprecated ``sunpy.io.cdf`` submodule, which is not intended to be user facing. (`#7240 <https://github.com/sunpy/sunpy/pull/7240>`__)
- Removed the deprecated ``sunpy.io.jp2``, which is not intended to be user facing. (`#7241 <https://github.com/sunpy/sunpy/pull/7241>`__)
- Removed the deprecated ``sunpy.io.file_tools``, which is not intended to be user facing. (`#7242 <https://github.com/sunpy/sunpy/pull/7242>`__)
- The deprecated ``sunpy.data.download_sample_data()`` has been removed
  Use :func:`sunpy.data.sample.download_all` instead. (`#7250 <https://github.com/sunpy/sunpy/pull/7250>`__)

Internal Changes
----------------

- Removed the Binder configuration and link in README.
  This is because the configuration was untested, and does not currently work. (`#7062 <https://github.com/sunpy/sunpy/pull/7062>`__)
- Add a Dependabot config file to auto-update GitHub action versions. (`#7068 <https://github.com/sunpy/sunpy/pull/7068>`__)
- Add tests to check whether various `~sunpy.map.Map` methods preserve laziness when operating on Maps backed by a `dask.array.Array`. (`#7100 <https://github.com/sunpy/sunpy/pull/7100>`__)
- Added missing support to find GOES-18 XRS data in `~sunpy.net.dataretriever.XRSClient`. (`#7108 <https://github.com/sunpy/sunpy/pull/7108>`__)
- Raise an error with a helpful message when :meth:`sunpy.map.GenericMap.plot` is called with a non-boolean value for the ``annotate`` keyword, because the user is probably trying to specify the axes. (`#7163 <https://github.com/sunpy/sunpy/pull/7163>`__)
- Fixed our ASDF manifest having the incorrect ID. (`#7282 <https://github.com/sunpy/sunpy/pull/7282>`__)
- Fix example formatting in a few asdf schemas. (`#7292 <https://github.com/sunpy/sunpy/pull/7292>`__)
- Pinned the ``drms`` requirement to ``< 0.7`` to avoid breaking changes in ``drms`` version 0.7. (`#7308 <https://github.com/sunpy/sunpy/pull/7308>`__)


5.0.0 (2023-06-14)
==================

Breaking Changes
----------------

- `~sunpy.net.dataretriever.XRSClient` now provides the re-processed GOES-XRS 8-15 data from NOAA.
  These files are now all NetCDF and not FITS files. (`#6737 <https://github.com/sunpy/sunpy/pull/6737>`__)
- Changed the output of :func:`sunpy.map.sample_at_coords` to return the sampled values as `~astropy.units.Quantity` with the appropriate units instead of merely numbers. (`#6882 <https://github.com/sunpy/sunpy/pull/6882>`__)


Deprecations
------------

- Using ``sunpy.map.header_helper.meta_keywords`` is deprecated.
  Please see :ref:`Meta Keywords Table` for the list of metadata keywords used by `~sunpy.map.Map`. (`#6743 <https://github.com/sunpy/sunpy/pull/6743>`__)
- The utility function ``sunpy.map.extract_along_coord`` is deprecated.
  Use :func:`sunpy.map.pixelate_coord_path`, and then pass its output to :func:`sunpy.map.sample_at_coords`. (`#6840 <https://github.com/sunpy/sunpy/pull/6840>`__)
- Parsing SDO/EVE level 0CS average files is deprecated, and will be removed in sunpy 6.0.
  Parsing this data is untested, and we cannot find a file to test it with.
  If you know where level 0CS 'averages' files can be found, please get in touch at https://community.openastronomy.org/c/sunpy/5. (`#6857 <https://github.com/sunpy/sunpy/pull/6857>`__)
- Fully deprecated ``sunpy.database``, with an expected removal version of sunpy 6.0. (`#6869 <https://github.com/sunpy/sunpy/pull/6869>`__)
- ``sunpy.io.cdf``, ``sunpy.io.file_tools`` and ``sunpy.io.jp2`` sub-modules have been deprecated, and will be removed in version 5.1.
  This because they are designed for internal use only, and removing it from the public API gives the developers more flexibility to modify it without impacting users. (`#6895 <https://github.com/sunpy/sunpy/pull/6895>`__)


New Features
------------

- A pure Python ``sunpy`` wheel is now published on PyPI with each release.
  ``pip`` will now default to installing the pure Python wheel instead of the source distribution on platforms other than Linux (x86-64) and macOS (x86-64 and ARM64).
  This should mean simpler and faster installs on such platforms, which includes the Raspberry Pi as well as some cloud computing services.

  This wheel does not contain the ``sunpy.io.ana`` compiled extension.
  If you need this extension (not available on Windows) you can install the ``sunpy`` source distribution with ``pip install --no-binary sunpy "sunpy[all]"``. (`#6175 <https://github.com/sunpy/sunpy/pull/6175>`__)
- Added three tutorials which replicate `~sunpy.map.CompositeMap` functionality (:ref:`sphx_glr_generated_gallery_plotting_AIA_HMI_composite.py`, :ref:`sphx_glr_generated_gallery_plotting_masked_composite_plot.py`, :ref:`sphx_glr_generated_gallery_plotting_three_map_composite.py`). (`#6459 <https://github.com/sunpy/sunpy/pull/6459>`__)
- `~sunpy.map.GenericMap.exposure_time` now looks for the exposure time in the ``XPOSURE`` key first
  and then the ``EXPTIME`` key. (`#6557 <https://github.com/sunpy/sunpy/pull/6557>`__)
- `~sunpy.map.header_helper.make_fitswcs_header` now includes the keyword argument ``detector`` for setting the
  ``DETECTOR`` FITS keyword in the resulting header. (`#6558 <https://github.com/sunpy/sunpy/pull/6558>`__)
- Adds two tutorials that demonstrate how to use LASCO data in overlaying maps (:ref:`sphx_glr_generated_gallery_plotting_lasco_overlay.py`) and how to create a custom mask for a LASCO C2 image (:ref:`sphx_glr_generated_gallery_map_lasco_mask.py`). (`#6576 <https://github.com/sunpy/sunpy/pull/6576>`__)
- Able to run the ``sunpy`` tests doing ``python -m sunpy.tests.self_test``. (`#6600 <https://github.com/sunpy/sunpy/pull/6600>`__)
- Able to detect gzip-compressed FITS files even if they don't have the ``.gz`` extension in the filename.
  ``sunpy.io.detect_filetype`` now looks for the right file signature while checking
  for gzipped FITS files. (`#6693 <https://github.com/sunpy/sunpy/pull/6693>`__)
- Added ``AttrAnd`` and ``AttrOr`` to the namespace in ``sunpy.net.attrs``.
  This allows users to to avoid ``|`` or ``&`` when creating a query a larger query. (`#6708 <https://github.com/sunpy/sunpy/pull/6708>`__)
- `~sunpy.net.dataretriever.SUVIClient` now provides GOES-18 SUVI data. (`#6737 <https://github.com/sunpy/sunpy/pull/6737>`__)
- The minimum required versions of several core dependencies have been updated:

  - Python 3.9
  - astropy 5.0.1
  - numpy 1.21.0

  The minimum required versions of these optional dependencies has also been updated:

  - Matplotlib 3.5.0
  - dask 2021.4.0
  - pandas 1.2.0
  - scikit-image 0.18.0
  - scipy 1.7.0 (`#6742 <https://github.com/sunpy/sunpy/pull/6742>`__)
- Added the utility function :func:`sunpy.map.pixelate_coord_path` to fully pixelate a coordinate path according to the pixels of a given map. (`#6840 <https://github.com/sunpy/sunpy/pull/6840>`__)
- The minimum version of h5netcdf required by sunpy has been bumped to version 0.11.0. (`#6859 <https://github.com/sunpy/sunpy/pull/6859>`__)
- Able to download files from REST/TAP Data Providers from the VSO. (`#6887 <https://github.com/sunpy/sunpy/pull/6887>`__)
- Adding data unit into html repr for `sunpy.map.Map` (`#6902 <https://github.com/sunpy/sunpy/pull/6902>`__)
- Joined ``HISTORY`` keys with newline characters when parsing ``HISTORY`` cards from
  FITS header. (`#6911 <https://github.com/sunpy/sunpy/pull/6911>`__)
- Added the ability to query for the GOES-XRS 1 minute average data with the `.XRSClient`. (`#6925 <https://github.com/sunpy/sunpy/pull/6925>`__)
- Increased minimum version of `parfive` to 2.0.0.

  We are aware the change in the ``parfive`` minimum version is a release earlier than our dependency policy allows for.
  However, due to significant issues that ``parfive`` v2.0.0 solves and changes to remote servers, we have decided to increase it to improve the user experience when downloading files. (`#6942 <https://github.com/sunpy/sunpy/pull/6942>`__)


Bug Fixes
---------

- Fixed the incorrect calculation in :func:`~sunpy.map.header_helper.make_fitswcs_header` of the rotation matrix from a rotation angle when the pixels are non-square. (`#6597 <https://github.com/sunpy/sunpy/pull/6597>`__)
- Return code from ``self_test`` is now non-zero if it stops due to missing dependencies. (`#6600 <https://github.com/sunpy/sunpy/pull/6600>`__)
- Fixed an issue with loading old EIT fits files with `sunpy.map.Map` where the date could not be parsed. (`#6605 <https://github.com/sunpy/sunpy/pull/6605>`__)
- Fixed a bug where the `~sunpy.map.GenericMap.exposure_time` returned ``None`` when the exposure
  time key was set to zero. (`#6637 <https://github.com/sunpy/sunpy/pull/6637>`__)
- Fixed a bug that prevented specifying a `~astropy.coordinates.BaseCoordinateFrame` (as opposed to a `~astropy.coordinates.SkyCoord`) to :meth:`sunpy.map.GenericMap.draw_quadrangle`. (`#6648 <https://github.com/sunpy/sunpy/pull/6648>`__)
- HMI JPEG2000 files from Helioviewer could not be loaded due to a bug in setting the plotting normalization.
  This has been fixed. (`#6710 <https://github.com/sunpy/sunpy/pull/6710>`__)
- The ``data_manager`` was not raising failed downloads correctly and would continue as if the file existed locally.
  Now it will raise any errors from ``parfive``. (`#6711 <https://github.com/sunpy/sunpy/pull/6711>`__)
- `~sunpy.map.sources.XRTMap` will now set the unit for XRT files if the ``BUNIT`` key is missing. (`#6725 <https://github.com/sunpy/sunpy/pull/6725>`__)
- `~sunpy.net.dataretriever.XRSClient` update use the new url for which the GOES-XRS 8-15 data is provided by NOAA. (`#6737 <https://github.com/sunpy/sunpy/pull/6737>`__)
- Updated ``sunpy.database`` to be compatible with ``SQLAlchemy`` versions >=2.0 (`#6749 <https://github.com/sunpy/sunpy/pull/6749>`__)
- When using ``autoalign=True`` when plotting maps, the result was misaligned by half a pixel. (`#6796 <https://github.com/sunpy/sunpy/pull/6796>`__)
- :meth:`sunpy.map.GenericMap.submap` can now handle a `~astropy.coordinates.BaseCoordinateFrame` as input. (`#6820 <https://github.com/sunpy/sunpy/pull/6820>`__)
- Multi-line ``HISTORY`` and ``COMMENT`` keys metadata dictionaries are now correctly split into
  multiple history and comment cards when writing a FITS file. (`#6911 <https://github.com/sunpy/sunpy/pull/6911>`__)
- Pass in "max_splits" to Parfive to prevent multi connections to JSOC for JSOC only queries. (`#6921 <https://github.com/sunpy/sunpy/pull/6921>`__)
- When converting an `astropy.wcs.WCS` object to a solar coordinate frame the
  ``DATE-AVG`` key will be used before the ``DATE-OBS`` key, previously only
  ``DATE-OBS`` was checked. (`#6995 <https://github.com/sunpy/sunpy/pull/6995>`__)
- `sunpy.map.GenericMap.rotation_matrix` now applies the default values if any FITS rotation matrix keywords are missing from the header. (`#7004 <https://github.com/sunpy/sunpy/pull/7004>`__)
- Modified :func:`sunpy.io.special.srs.read_srs` to correctly handle uppercase SRS files and supplementary sections occurring after the main data sections (I, IA, II). (`#7035 <https://github.com/sunpy/sunpy/pull/7035>`__)


Documentation
-------------

- Added an example of how to search for multiple wavelengths attributes for AIA data using `sunpy.net.attrs.AttrOr`. (`#6501 <https://github.com/sunpy/sunpy/pull/6501>`__)
- Added `sunpy.map.PixelPair` to the reference documentation. (`#6620 <https://github.com/sunpy/sunpy/pull/6620>`__)
- Split the installation docs into a new Installation tutorial, and an installation guide. (`#6639 <https://github.com/sunpy/sunpy/pull/6639>`__)
- Added an example (:ref:`sphx_glr_generated_gallery_time_series_goes_xrs_nrt_data.py`) to download GOES NRT data and load it into `~sunpy.timeseries.TimeSeries`. (`#6744 <https://github.com/sunpy/sunpy/pull/6744>`__)
- Added an example gallery (:ref:`sphx_glr_generated_gallery_acquiring_data_querying_and_loading_SHARP_data.py`) for querying SHARP data and loading it into a `~sunpy.map.Map`. (`#6757 <https://github.com/sunpy/sunpy/pull/6757>`__)
- Added an example (:ref:`sphx_glr_generated_gallery_units_and_coordinates_ParkerSolarProbe_trajectory.py`) to plot the trajectory of Parker Solar Probe. (`#6771 <https://github.com/sunpy/sunpy/pull/6771>`__)
- Created a "Showcase" section of the gallery, which includes a new example (:ref:`sphx_glr_generated_gallery_showcase_where_is_stereo.py`) and a relocated example (:ref:`sphx_glr_generated_gallery_showcase_hmi_cutout.py`). (`#6781 <https://github.com/sunpy/sunpy/pull/6781>`__)
- Updated examples in the gallery to always explicitly create an Axes and use that for plotting, instead of using the Matplotlib pyplot API. (`#6822 <https://github.com/sunpy/sunpy/pull/6822>`__)
- Added an example (:ref:`sphx_glr_generated_gallery_map_masking_hmi.py`) of how to mask a HMI map based on the intensity of AIA. (`#6825 <https://github.com/sunpy/sunpy/pull/6825>`__)
- Added an example to blend two maps using ``mplcairo``. (`#6835 <https://github.com/sunpy/sunpy/pull/6835>`__)
- Changed the reprojecting images to different observers example (:ref:`sphx_glr_generated_gallery_map_transformations_reprojection_different_observers.py`) to avoid using custom wcs headers where possible. (`#6853 <https://github.com/sunpy/sunpy/pull/6853>`__)
- Added a note in examples :ref:`sphx_glr_generated_gallery_map_transformations_autoalign_aia_hmi.py` and :ref:`sphx_glr_generated_gallery_map_transformations_reprojection_align_aia_hmi.py` suggesting to use :meth:`~sunpy.coordinates.Helioprojective.assume_spherical_screen` to retain off-disk HMI data. (`#6855 <https://github.com/sunpy/sunpy/pull/6855>`__)
- Moved the Helioviewer migration guide from the tutorial to guide section of the docs. (`#6868 <https://github.com/sunpy/sunpy/pull/6868>`__)
- Moved the plotting section of the tutorial into the map section of the tutorial. (`#6870 <https://github.com/sunpy/sunpy/pull/6870>`__)
- Reorganized "Units" section of the Tutorial into smaller sections and added a section about
  unit equivalencies. (`#6879 <https://github.com/sunpy/sunpy/pull/6879>`__)
- Added clarifying detail (in the `~sunpy.time.TimeUTime` docstring) for how the ``utime`` time format handles seconds on a day with a leap second. (`#6894 <https://github.com/sunpy/sunpy/pull/6894>`__)
- Fixed a series of broken URLS and typos in examples and documentation strings. (`#6903 <https://github.com/sunpy/sunpy/pull/6903>`__)
- Improved the time tutorial. (`#6920 <https://github.com/sunpy/sunpy/pull/6920>`__)
- Add a "how-to" guide section to the documentation. (`#6926 <https://github.com/sunpy/sunpy/pull/6926>`__)
- Redesigned the landing page to highlight the different sections of the documentation. (`#6938 <https://github.com/sunpy/sunpy/pull/6938>`__)
- Significantly revised and improved the :ref:`sunpy-tutorial-maps` part of the tutorial.
  This included moving the section on custom maps to the :ref:`sunpy-how-to-index` section (see :ref:`sunpy-how-to-create-a-map`). (`#6944 <https://github.com/sunpy/sunpy/pull/6944>`__)
- Migrated example gallery entries for searching the VSO, using ``parse_time``, using the data manager, and using solar constants to the how-to guide. (`#6948 <https://github.com/sunpy/sunpy/pull/6948>`__)
- Reorganized some parts of the coordinates topic guide into multiple how-to guides. (`#6954 <https://github.com/sunpy/sunpy/pull/6954>`__)
- Move examples of how to create a Map from reference pages to a how-to guide. (`#6977 <https://github.com/sunpy/sunpy/pull/6977>`__)
- Cleaned up and simplified the :ref:`sunpy-tutorial-timeseries` section of the tutorial. (`#6990 <https://github.com/sunpy/sunpy/pull/6990>`__)
- Added a topic-guide to aid understanding the role, "rsun" plays in sunpy coordinate transformations and :meth:`sunpy.map.GenericMap.reproject_to`. (`#7000 <https://github.com/sunpy/sunpy/pull/7000>`__)
- Updated all of the sphinx anchors to be more consistent.
  This means that any use of the old anchors (intersphinx links to sunpy doc pages) will need to be updated. (`#7032 <https://github.com/sunpy/sunpy/pull/7032>`__)


Internal Changes
----------------

- When determining which VSO servers to use for queries, `.VSOClient` will now
  attempt to check if the cgi endpoint referenced by the WDSL file is accessible,
  and try the next endpoint if it can't be reached. This should mean that a small
  category of connection issues with the VSO are now automatically bypassed. (`#6362 <https://github.com/sunpy/sunpy/pull/6362>`__)


4.1.0 (2022-11-11)
==================

Breaking Changes
----------------

- Updated the sample data file, ``AIA_171_ROLL_IMAGE`` to be rice compressed instead of gzip compressed.
  This means that the data is now stored in the second HDU. (`#6221 <https://github.com/sunpy/sunpy/pull/6221>`__)


Deprecations
------------

- Passing positional arguments to all ``timeseries`` ``peek()`` methods
  is now deprecated, and will raise an error in sunpy 5.1. Pass the arguments
  with keywords (e.g. ``title='my plot title'``) instead. (`#6310 <https://github.com/sunpy/sunpy/pull/6310>`__)
- Using ``sunpy.timeseries.GenericTimeSeries.index``` is deprecated.
  Use `~sunpy.timeseries.GenericTimeSeries.time` to get an astropy Time object,
  or ``ts.to_dataframe().index`` to get the times as a pandas ``DataTimeIndex``. (`#6327 <https://github.com/sunpy/sunpy/pull/6327>`__)
- Deprecated the ``sunpy.visualization.limb`` module.
  The ``sunpy.visualization.limb.draw_limb`` function has been moved into
  `~sunpy.visualization.drawing` as :func:`~sunpy.visualization.drawing.limb`. (`#6332 <https://github.com/sunpy/sunpy/pull/6332>`__)
- The ``sunpy.net.helioviewer`` module is deprecated and will be removed in version 5.1.
  The Helioviewer Project now maintains a replacement Python library called `hvpy <https://hvpy.readthedocs.io/en/latest/>`__.
  As such, in consultation with the Helioviewer Project, we have decided to deprecate the ``HelioviewerClient`` class. (`#6404 <https://github.com/sunpy/sunpy/pull/6404>`__)
- Passing the ``algorithm``, ``return_footprint`` arguments as positional arguments is deprecated. Pass them as keyword arguments (e.g. ``..., return_footprint=True, ...``) instead. (`#6406 <https://github.com/sunpy/sunpy/pull/6406>`__)
- ``sunpy.data.download_sample_data()`` is now deprecated.
  Use :func:`sunpy.data.sample.download_all` instead. (`#6426 <https://github.com/sunpy/sunpy/pull/6426>`__)
- The sunpy.database module is no longer actively maintained and has a number of outstanding issues.
  It is anticipated that sunpy.database will be formally deprecated in sunpy 5.0 and removed in sunpy 6.0.
  If you are using sunpy.database and would like to see a replacement, please join the discussion thread at https://community.openastronomy.org/t/deprecating-sunpy-database/495. (`#6498 <https://github.com/sunpy/sunpy/pull/6498>`__)


Removals
--------

- The ``sunpy.io.fits`` sub-module has been removed, as it was designed for internal use.
  Use the `astropy.io.fits` module instead for more generic functionality to read FITS files. (`#6432 <https://github.com/sunpy/sunpy/pull/6432>`__)
- The ``sunpy.physics.solar_rotation`` sub-module has been removed, having been moved to `sunkit_image.coalignment`. (`#6433 <https://github.com/sunpy/sunpy/pull/6433>`__)
- Most of the `sunpy.visualization.animator` subpackage has been removed, with the exception of `~sunpy.visualization.animator.MapSequenceAnimator`
  It has been moved into the standalone `mpl-animators <https://pypi.org/project/mpl-animators>`_ package
  Please update your imports to replace ``sunpy.visualization.animator`` with ``mpl_animators``. (`#6434 <https://github.com/sunpy/sunpy/pull/6434>`__)
- Remove ``GenericMap.shift`` method and the ``GenericMap.shifted_value``.
  Use `~sunpy.map.GenericMap.shift_reference_coord` instead. (`#6437 <https://github.com/sunpy/sunpy/pull/6437>`__)
- ``sunpy.util.scraper`` has been removed. Use `sunpy.net.scraper` instead. (`#6438 <https://github.com/sunpy/sunpy/pull/6438>`__)
- ``sunpy.image.coalignment`` has been removed. Use `sunkit_image.coalignment` instead, which contains all the same functionality. (`#6440 <https://github.com/sunpy/sunpy/pull/6440>`__)
- :meth:`sunpy.map.GenericMap.draw_limb` can no longer be used to draw the limb on a non-WCS Axes plot. (`#6533 <https://github.com/sunpy/sunpy/pull/6533>`__)
- :meth:`sunpy.image.resample` no longer accepts "neighbour" as an interpolation method.
  Use "nearest" instead. (`#6537 <https://github.com/sunpy/sunpy/pull/6537>`__)
- :meth:`sunpy.image.transform.affine_transform` and :func:`sunpy.map.GenericMap.rotate` no longer accepts the ``use_scipy`` keyword. (`#6538 <https://github.com/sunpy/sunpy/pull/6538>`__)


New Features
------------

- Updated and expanded the HTML representation for `~sunpy.timeseries.TimeSeries`. (`#5951 <https://github.com/sunpy/sunpy/pull/5951>`__)
- When reading CDF files, any columns with a floating point data type now have their masked values converted to NaN. (`#5956 <https://github.com/sunpy/sunpy/pull/5956>`__)
- Add support for saving `~sunpy.map.GenericMap` as JPEG 2000 files. (`#6153 <https://github.com/sunpy/sunpy/pull/6153>`__)
- Add a function ``sunpy.map.extract_along_coord`` that, for a given set of coordinates,
  finds each array index that crosses the line traced by those coordinates and returns the value of the data
  array of a given map at those array indices. (`#6189 <https://github.com/sunpy/sunpy/pull/6189>`__)
- Three new maps have been added to the sample data from STEREO A and STEREO B at
  195 Angstrom, and AIA at 193 Angstrom. These images are from a time when
  the three spacecraft were equally spaced around the Sun, and therefore form
  near complete instantaneous coverage of the solar surface.

  Users upgrading to this version will find this three files download when they
  use the sample data for the first time. (`#6197 <https://github.com/sunpy/sunpy/pull/6197>`__)
- Added a SDO/AIA 1600 file of the Venus transit to the sunpy sample data. (`#6242 <https://github.com/sunpy/sunpy/pull/6242>`__)
- Created the `sunpy.visualization.drawing` module which includes
  new :func:`~sunpy.visualization.drawing.equator` and
  :func:`~sunpy.visualization.drawing.prime_meridian` functions. (`#6251 <https://github.com/sunpy/sunpy/pull/6251>`__)
- Expose GOES quality flags in order to allow filtering corrupt values when using the `~sunpy.timeseries.sources.goes.XRSTimeSeries`. (`#6260 <https://github.com/sunpy/sunpy/pull/6260>`__)
- All TimeSeries plotting methods now consistently set the same
  formatter and locator for the x-axis. (`#6264 <https://github.com/sunpy/sunpy/pull/6264>`__)
- :meth:`sunpy.timeseries.GenericTimeSeries.peek` now takes a ``title`` argument
  to set the title of the plot. (`#6304 <https://github.com/sunpy/sunpy/pull/6304>`__)
- Added the `sunpy.timeseries.GenericTimeSeries.time` property to get the times
  of a timeseries as a `~astropy.time.Time` object. (`#6327 <https://github.com/sunpy/sunpy/pull/6327>`__)
- Added the :ref:`sphx_glr_generated_gallery_plotting_plot_equator_prime_meridian.py` example to the Example Gallery. (`#6332 <https://github.com/sunpy/sunpy/pull/6332>`__)
- Added a new function :func:`sunpy.map.header_helper.make_heliographic_header` to help with generating FITS-WCS headers in Carrington or Stonyhurst coordinate systems that span the entire solar surface. (`#6415 <https://github.com/sunpy/sunpy/pull/6415>`__)
- Sample data files provided through `sunpy.data.sample` are now downloaded individually on demand rather than being all downloaded upon import of that module.
  To download all sample data files, call :func:`sunpy.data.sample.download_all`. (`#6426 <https://github.com/sunpy/sunpy/pull/6426>`__)
- `~.XRSTimeSeries` is now able to parse the primary detector information from the GOES-R XRS data if available. (`#6454 <https://github.com/sunpy/sunpy/pull/6454>`__)
- `sunpy.net.Scraper` now includes treats files as spanning a full interval equal to the smallest increment specified in the file pattern.
  For example, a pattern like ``"%Y.txt"`` that only contains a year specifier will be considered to span that full year.

  This means searches that fall entirely within the whole interval spanned by a pattern will return that file, where previously they did not.
  As an example, matching ``"%Y.txt"`` with ``TimeRange('2022-02-01', '2022-04-01')`` will now return ``["2022.txt"]`` where previously no files were returned. (`#6472 <https://github.com/sunpy/sunpy/pull/6472>`__)
- Implemented site configuration for sunpyrc, and modified documentation for sunpy customization. (`#6478 <https://github.com/sunpy/sunpy/pull/6478>`__)
- :func:`~sunpy.map.header_helper.make_fitswcs_header` now includes the keyword argument ``unit`` for setting the
  ``BUNIT`` FITS keyword in the resulting header.
  This will take precedence over any unit information attached to ``data``. (`#6499 <https://github.com/sunpy/sunpy/pull/6499>`__)
- If the ``data`` argument to :func:`~sunpy.map.header_helper.make_fitswcs_header` is an `~astropy.units.Quantity`,
  the associated unit will be used to set the ``BUNIT`` FITS keyword in the resulting header. (`#6499 <https://github.com/sunpy/sunpy/pull/6499>`__)
- Added a 304 sample data file called ``AIA_304_IMAGE``. (`#6546 <https://github.com/sunpy/sunpy/pull/6546>`__)


Bug Fixes
---------

- Fix a bug that prevented EUI maps with missing wavelength metadata loading. (`#6199 <https://github.com/sunpy/sunpy/pull/6199>`__)
- The `sunpy.net.dataretriever.sources.noaa.SRSClient` was not correctly setting the passive mode for FTP connection resulting in a permission error.
  This has been fixed. (`#6256 <https://github.com/sunpy/sunpy/pull/6256>`__)
- Fixed `~sunpy.timeseries.sources.XRSTimeSeries` inability to read leap-second files for GOES.
  It floors the leap-second timestamp to be ``59.999``, so that Python datetime does not raise an exception. (`#6262 <https://github.com/sunpy/sunpy/pull/6262>`__)
- Changed the default scaling for `~sunpy.map.sources.EUIMap` from a linear stretch to a asinh stretch.

  To revert to the previous linear stretch do the following::

       from astropy.visualization import ImageNormalize, LinearStretch
       euimap.plot_settings["norm"] = ImageNormalize(stretch=LinearStretch()) (`#6285 <https://github.com/sunpy/sunpy/pull/6285>`__)
- Fixed bugs when working with a coordinate frame where the observer is specified in `~sunpy.coordinates.frames.HeliographicStonyhurst` with a Cartesian representation, which is equivalent to Heliocentric Earth Equatorial (HEEQ).
  Now, the observer will always be converted to spherical representation when the coordinate frame is created. (`#6311 <https://github.com/sunpy/sunpy/pull/6311>`__)
- Fixed an error when Fido returns zero results from the VSO
  and some results from at least one other data source. This
  (now fixed) error is only present when using numpy version >= 1.23. (`#6318 <https://github.com/sunpy/sunpy/pull/6318>`__)
- If a level 1 XRT file does not specify the heliographic longitude of the spacecraft,
  a silent assumption is made that the spacecraft is at zero Stonyhurst
  heliographic longitude (i.e., the same longitude as Earth). (`#6333 <https://github.com/sunpy/sunpy/pull/6333>`__)
- The sample data retry was failing under parfive 2.0.0. (`#6334 <https://github.com/sunpy/sunpy/pull/6334>`__)
- Fixed bug that prevented `~sunpy.coordinates.metaframes.RotatedSunFrame` instances from being pickled. (`#6342 <https://github.com/sunpy/sunpy/pull/6342>`__)
- Fix a bug in loading `.XRSTimeSeries` due to unsupported quality flag column names. (`#6410 <https://github.com/sunpy/sunpy/pull/6410>`__)
- Adds units (dimensionless units) to the quality columns in `.XRSTimeSeries`. (`#6423 <https://github.com/sunpy/sunpy/pull/6423>`__)
- Refactored `~sunpy.map.sources.SXTMap` to use ITRS observer coordinate information
  in header rather than incorrect HGS keywords.
  The `~sunpy.map.sources.SXTMap` also now uses the default ``dsun`` property as this
  information can be derived from the (now corrected) observer coordinate. (`#6436 <https://github.com/sunpy/sunpy/pull/6436>`__)
- In `sunpy.map.GenericMap.coordinate_system` and `sunpy.map.GenericMap.date`, the default values
  will now be used if the expected key(s) used to derive those properties are empty.
  Previously, empty values of these keys were not treated as missing and thus the default values
  were not correctly filled in. (`#6436 <https://github.com/sunpy/sunpy/pull/6436>`__)
- Fixed a bug where the observer coordinate was incorrectly determined for `~sunpy.map.sources.KCorMap`. (`#6447 <https://github.com/sunpy/sunpy/pull/6447>`__)
- Trying to download an empty search response from the JSOC now results in an empty results object.
  Previously the results object contained the path to the sunpy download directory. (`#6449 <https://github.com/sunpy/sunpy/pull/6449>`__)
- Removed an error when searching CDAWEB using `sunpy.net.Fido` and no results are returned.
  An empty response table is now returned. (`#6450 <https://github.com/sunpy/sunpy/pull/6450>`__)
- Fix a bug to parse the GOES "observatory" number in `~.XRSTimeSeries` for GOES 13, 14, 15 and for the 1 minute GOES-R data. (`#6451 <https://github.com/sunpy/sunpy/pull/6451>`__)
- Changed the default scaling for `~sunpy.map.sources.XRTMap` from a linear stretch to `~astropy.visualization.LogStretch`.

  To revert to the previous linear stretch do the following::

       from astropy.visualization import ImageNormalize, LinearStretch
       xrtmap.plot_settings["norm"] = ImageNormalize(stretch=LinearStretch()) (`#6480 <https://github.com/sunpy/sunpy/pull/6480>`__)
- Fix the ``detector`` property of `~sunpy.map.sources.SOTMap` to return "SOT". (`#6480 <https://github.com/sunpy/sunpy/pull/6480>`__)
- The right-hand y-axis of the GOES-XRS timeseries plots with labelled flare classes
  now automatically scales with the left-hand y-axis. (`#6486 <https://github.com/sunpy/sunpy/pull/6486>`__)
- Add support for Python 3.11.

  The deprecated "cgi.parse_header" is now available as
  `sunpy.util.net.parse_header`. (`#6512 <https://github.com/sunpy/sunpy/pull/6512>`__)
- Fixed the metadata handling of :meth:`~sunpy.map.GenericMap.resample` and :meth:`~sunpy.map.GenericMap.superpixel` so that the CDELTi values are scaled and the PCi_j matrix (if used) is modified in the correct manner for asymmetric scaling.
  The previous approach of having the PCi_j matrix store all of the scaling resulted in non-intuitive behaviors when accessing the `~sunpy.map.GenericMap.scale` and `~sunpy.map.GenericMap.rotation_matrix` properties, and when de-rotating a map via :meth:`~sunpy.map.GenericMap.rotate`. (`#6571 <https://github.com/sunpy/sunpy/pull/6571>`__)
- Fixd a bug with the `sunpy.map.GenericMap.scale` property for maps containing only the CDij matrix where the scale was not being determined from the CDij matrix. (`#6573 <https://github.com/sunpy/sunpy/pull/6573>`__)
- Fixed a bug with the `sunpy.map.GenericMap.rotation_matrix` property for maps using the CDij matrix formulism where the rotation matrix would be calculated incorrectly for non-square pixels. (`#6573 <https://github.com/sunpy/sunpy/pull/6573>`__)
- Fixed a bug where :func:`~sunpy.time.parse_time` would always disregard the remainder of a time string starting with the final period if it was followed by only zeros, which could affect the parsing of the time string. (`#6581 <https://github.com/sunpy/sunpy/pull/6581>`__)


Documentation
-------------

- Improved annotations in the SRS active regions plotting example. (`#6196 <https://github.com/sunpy/sunpy/pull/6196>`__)
- Updated gallery examples that use STEREO data to use sample data instead
  of searching for and downloading data via Fido. (`#6197 <https://github.com/sunpy/sunpy/pull/6197>`__)
- Added the current bugfix release policy to the docs. (`#6336 <https://github.com/sunpy/sunpy/pull/6336>`__)
- The :ref:`sunpy-tutorial-maps` and :ref:`sunpy-tutorial-timeseries` have been reviewed and updated. (`#6345 <https://github.com/sunpy/sunpy/pull/6345>`__)
- Adds a pull request check list to the Developer's Guide. (`#6346 <https://github.com/sunpy/sunpy/pull/6346>`__)
- Improved the plotting guide. (`#6430 <https://github.com/sunpy/sunpy/pull/6430>`__)
- Slight improvements to the downloading data with Fido part of the guide. (`#6444 <https://github.com/sunpy/sunpy/pull/6444>`__)
- Split the units and coordinate guides on to separate pages, and made minor improvements to them. (`#6462 <https://github.com/sunpy/sunpy/pull/6462>`__)
- Added a how-to guide ``conda_for_dependencies`` for using ``conda`` to set up an environment with the complete set of dependencies to use all optional features, build the documentation, and/or run the full test suite.
  The guide also describes how best to have an editable installation of ``sunpy`` in this environment. (`#6524 <https://github.com/sunpy/sunpy/pull/6524>`__)


Internal Changes
----------------

- Added a ``columns`` keyword to each plot method for all `sunpy.timeseries.GenericTimeSeries` sources. (`#6056 <https://github.com/sunpy/sunpy/pull/6056>`__)
- Added a script in the ``sunpy/tools`` that will update all the Python libraries in ``sunpy/extern``. (`#6127 <https://github.com/sunpy/sunpy/pull/6127>`__)
- Added automatic conversion of unit strings in CDF files to astropy unit objects for the following instruments: PSP/ISOIS, SOHO/CELIAS, SOHO/COSTEP-EPHIN, and SOHO/ERNE. (`#6159 <https://github.com/sunpy/sunpy/pull/6159>`__)
- Add an environment variable ``SUNPY_NO_BUILD_ANA_EXTENSION`` which when present
  will cause sunpy to not compile the ANA C extension when building from source. (`#6166 <https://github.com/sunpy/sunpy/pull/6166>`__)
- ``sunpy`` now uses the `Limited Python API <https://docs.python.org/3/c-api/stable.html>`__.
  Therefore, one binary distribution (wheel) per platform is now published and it is compatible with all Python versions ``sunpy`` supports. (`#6171 <https://github.com/sunpy/sunpy/pull/6171>`__)
- Add support for upcoming parfive 2.0 release. (`#6243 <https://github.com/sunpy/sunpy/pull/6243>`__)
- The primary sample-data URL will be changing from ``https://github.com/sunpy/sample-data/raw/master/sunpy/v1/`` to ``https://github.com/sunpy/data/raw/main/sunpy/v1/``.
  We expect GitHub to redirect from the old URL for sometime but will eventually expire it.
  The ``data.sunpy.org`` mirror will continue to be available. (`#6289 <https://github.com/sunpy/sunpy/pull/6289>`__)
- Add support for downloading sample data from more than two mirror locations. (`#6295 <https://github.com/sunpy/sunpy/pull/6295>`__)
- Timeseries data sources can now set the ``_peek_title`` class attribute
  to set the default plot title produced when ``.peek()`` is called and the user
  does not provide a custom title. (`#6304 <https://github.com/sunpy/sunpy/pull/6304>`__)
- All internal code for limb drawing now uses :func:`~sunpy.visualization.drawing.limb`. (`#6332 <https://github.com/sunpy/sunpy/pull/6332>`__)
- Add maintainer documentation on the backport bot (`#6355 <https://github.com/sunpy/sunpy/pull/6355>`__)
- Switched to using the standard matrix-multiplication operator (available in Python 3.5+) instead of a custom function. (`#6376 <https://github.com/sunpy/sunpy/pull/6376>`__)
- Fixed a colormap deprecation warning when importing the sunpy colormaps
  with Matplotlib 3.6. (`#6379 <https://github.com/sunpy/sunpy/pull/6379>`__)
- Removed custom tick label rotation from Lyra, EVE, and Norh timeseries sources, and grid drawing from NOAA and RHESSI sources. (`#6385 <https://github.com/sunpy/sunpy/pull/6385>`__)
- Added tests and test data for `~sunpy.map.sources.SXTMap` (`#6436 <https://github.com/sunpy/sunpy/pull/6436>`__)
- Fixed a bug where the private attribute ``_default_observer_coordinate`` for `~sunpy.map.GenericMap` was being used even when there was sufficient observer metadata in the header. (`#6447 <https://github.com/sunpy/sunpy/pull/6447>`__)
- Tidy the GOES XRSTimesSeries tests and add two new XRS files to test. (`#6460 <https://github.com/sunpy/sunpy/pull/6460>`__)
- Added a pre-commit hook for `codespell
  <https://github.com/codespell-project/codespell>`__, and applied
  spelling fixes throughout the package. (`#6574 <https://github.com/sunpy/sunpy/pull/6574>`__)


v4.0.0 (2022-05-06)
===================

Breaking Changes
----------------

- When rotating images using the SciPy rotation method, the default behavior is now to clip the output range to the input range, which matches the default behavior of the scikit-image rotation method. (`#5867 <https://github.com/sunpy/sunpy/pull/5867>`__)
- Any NaNs are now preserved by :func:`sunpy.image.transform.affine_transform` and :meth:`sunpy.map.GenericMap.rotate`. (`#5867 <https://github.com/sunpy/sunpy/pull/5867>`__)
- :func:`sunpy.image.transform.affine_transform` and :meth:`sunpy.map.GenericMap.rotate` now default to using SciPy for rotation instead of scikit-image, so rotation results may be slightly different. (`#5867 <https://github.com/sunpy/sunpy/pull/5867>`__)
- The math convenience methods of `sunpy.map.GenericMap` - :meth:`~sunpy.map.GenericMap.max`, :meth:`~sunpy.map.GenericMap.mean`, :meth:`~sunpy.map.GenericMap.min`, and , :meth:`~sunpy.map.GenericMap.std` - now ignore NaNs in the image data. (`#5867 <https://github.com/sunpy/sunpy/pull/5867>`__)
- :func:`sunpy.image.transform.affine_transform` and :meth:`sunpy.map.GenericMap.rotate` now default to using NaN instead of zero for the ``missing`` value, the value used for pixels in the output array that have no corresponding pixel in the input array.
  To obtain the previous behavior, ``missing`` should be explicitly specified as zero. (`#5867 <https://github.com/sunpy/sunpy/pull/5867>`__)
- The `.JSOCClient` and every `sunpy.net.dataretriever.GenericClient` was passing all ``**kwargs`` to `parfive.Downloader.enqueue_file`, this was unintended and has been removed. (`#6052 <https://github.com/sunpy/sunpy/pull/6052>`__)
- Changed the default interpolation order for :meth:`sunpy.map.GenericMap.rotate` from 4 to 3, with the precise meaning of these interpolation orders depending on the selected rotation method.
  For the default rotation method, which uses :func:`scipy.ndimage.affine_transform`, this changes the default interpolation from biquartic to bicubic, which reduces the computation time without reducing the quality of the output below what a typical user needs. (`#6089 <https://github.com/sunpy/sunpy/pull/6089>`__)


Deprecations
------------

- Deprecate ``sunpy.image.coalignment`` as the code has now been moved to
  `sunkit_image.coalignment` with an identical API.
  This module will be removed in sunpy 4.1. (`#5957 <https://github.com/sunpy/sunpy/pull/5957>`__)
- The ``sunpy.map.GenericMap.shift`` method has been renamed to
  `sunpy.map.GenericMap.shift_reference_coord` and
  ``shift`` has been deprecated. (`#5977 <https://github.com/sunpy/sunpy/pull/5977>`__)
- The ``sunpy.map.GenericMap.shifted_value`` property has been deprecated.
  Modifications to the reference coordinate can be found in the
  ``CRVAL1`` and ``CRVAL2`` keys of ``sunpy.map.GenericMap.meta.modified_items``. (`#5977 <https://github.com/sunpy/sunpy/pull/5977>`__)
- The ``sunpy.io.fits`` module is deprecated, as it was designed for internal use
  only. Use the `astropy.io.fits` module instead for more generic functionality
  to read FITS files. (`#5983 <https://github.com/sunpy/sunpy/pull/5983>`__)
- ``sunpy.physics.solar_rotation.mapsequence_solar_derotate`` is deprecated and will be removed in version 4.1.
  This function has been moved to ``sunkit_image.coalignment.mapsequence_coalign_by_rotation`` and has an identical API and functionality. (`#6031 <https://github.com/sunpy/sunpy/pull/6031>`__)
- ``sunpy.physics.solar_rotation.calculate_solar_rotate_shift`` is deprecated and will be removed in version 4.1.
  This function has been moved to ``sunkit_image.coalignment.calculate_solar_rotate_shift`` and has an identical API and functionality. (`#6031 <https://github.com/sunpy/sunpy/pull/6031>`__)
- Deprecated using `sunpy.map.GenericMap.draw_limb` on an Axes that is not a
  WCSAxes. (`#6079 <https://github.com/sunpy/sunpy/pull/6079>`__)


New Features
------------

- Added support for Python 3.10 (`#5568 <https://github.com/sunpy/sunpy/pull/5568>`__)
- Added support for ``"%Y.%m.%d_%H:%M:%S_UTC"`` and ``"%Y.%m.%d_%H:%M:%S"`` time formats in `sunpy.time.parse_time`. (`#5647 <https://github.com/sunpy/sunpy/pull/5647>`__)
- The ``rsun`` argument to :func:`~sunpy.map.header_helper.get_observer_meta` is now
  optional. (`#5655 <https://github.com/sunpy/sunpy/pull/5655>`__)
- Added the :meth:`~sunpy.net.base_client.QueryResponseTable.total_size`, which
  estimates the total size of the results from a Fido query. If this is supported
  by a client, the total size is printed alongside the results.

  To add support for this in external clients, make sure one column contains
  the individual filesizes as `~astropy.units.Quantity`, and set the
  ``size_column`` class attribute to the name of this column. (`#5659 <https://github.com/sunpy/sunpy/pull/5659>`__)
- Added the ability to specify the use of Carrington coordinates with
  :meth:`sunpy.map.GenericMap.draw_grid`. (`#5703 <https://github.com/sunpy/sunpy/pull/5703>`__)
- Printing a `.MetaDict`  will now show each entry on a new line. (`#5765 <https://github.com/sunpy/sunpy/pull/5765>`__)
- Removed support for Python 3.7. (`#5773 <https://github.com/sunpy/sunpy/pull/5773>`__)
- The 'event_endtime', 'event_starttime' and 'event_peaktime' columns in a HEK
  query are now returned as `~astropy.time.Time` objects. Previously they were
  timestamp strings. (`#5806 <https://github.com/sunpy/sunpy/pull/5806>`__)
- Added a helpful warning message when converting a 2D Helioprojective coordinate will return all NaNs. (`#5817 <https://github.com/sunpy/sunpy/pull/5817>`__)
- The colorbar limits on HMI magnetic field maps are now automatically
  set to be symmetric about zero. (`#5825 <https://github.com/sunpy/sunpy/pull/5825>`__)
- Added a ``clip`` keyword to :func:`sunpy.image.transform.affine_transform` and :meth:`sunpy.map.GenericMap.rotate` to enable or disable whether the range of the output image is clipped to the range of the input range. (`#5867 <https://github.com/sunpy/sunpy/pull/5867>`__)
- Created the decorator :func:`sunpy.image.transform.add_rotation_function` for registering new rotation functions for use by :func:`sunpy.image.transform.affine_transform` and :meth:`sunpy.map.GenericMap.rotate`. (`#5867 <https://github.com/sunpy/sunpy/pull/5867>`__)
- `sunpy.image.transform.affine_transform` and :meth:`sunpy.map.GenericMap.rotate`
  have both had their ``use_scipy`` arguments deprecated. Instead use the new
  ``method`` argument to select from the available rotation methods. (`#5916 <https://github.com/sunpy/sunpy/pull/5916>`__)
- Added a Maxwell unit and any places where a conversion to Gauss occurs has been removed. (`#5998 <https://github.com/sunpy/sunpy/pull/5998>`__)
- Add a basic HTML representation for `~sunpy.timeseries.TimeSeries`. (`#6032 <https://github.com/sunpy/sunpy/pull/6032>`__)
- The minimum supported asdf version has been increased to 2.8.0 to allow future
  compatibility with the breaking changes planned for asdf 3.0.
  In addition to this the `asdf-astropy <https://github.com/astropy/asdf-astropy>`__
  package is now required to serialise and deserialise the sunpy coordinate frame
  classes to ASDF. (`#6057 <https://github.com/sunpy/sunpy/pull/6057>`__)
- Added the option to rotate using `OpenCV <https://opencv.org>`__ when using :func:`sunpy.image.transform.affine_transform` or :meth:`sunpy.map.GenericMap.rotate` by specifying ``method='cv2'``.
  The OpenCV Python package must be installed on the system. (`#6089 <https://github.com/sunpy/sunpy/pull/6089>`__)


Bug Fixes
---------

- Fixed reading CDF files when a column has no entries. If this is the case the
  column will be ignored, and a message logged at DEBUG level. (`#5664 <https://github.com/sunpy/sunpy/pull/5664>`__)
- Fixed the units of `sunpy.map.sources.HMISynopticMap.scale` and
  `sunpy.map.sources.MDISynopticMap.scale`. (`#5682 <https://github.com/sunpy/sunpy/pull/5682>`__)
- Fixed a bug where custom values in the ``plot_settings`` dictionary were not being propagated
  to new map instances created when calling map methods (e.g. ``.submap``). (`#5687 <https://github.com/sunpy/sunpy/pull/5687>`__)
- Added automatic conversion of some common but non-standard unit strings in CDF
  files to astropy unit objects. If sunpy does not recognise the unit string for
  a particular column, units of ``u.dimensionless_unscaled`` are applied to that
  column and a warning raised.

  If you think a given unit should not be dimensionless and support should be
  added for it in sunpy, please raise an issue at
  https://github.com/sunpy/sunpy/issues. (`#5692 <https://github.com/sunpy/sunpy/pull/5692>`__)
- The default ``id_type`` in :func:`sunpy.coordinates.get_horizons_coord` is now
  `None` to match the default ``id_type`` in astroquery 0.4.4, which will search
  major bodies first, and if no major bodies are found, then search small bodies.
  For older versions of astroquery the default ``id_type`` used by
  :func:`~sunpy.coordinates.get_horizons_coord` is still ``'majorbody'``. (`#5707 <https://github.com/sunpy/sunpy/pull/5707>`__)
- In consultation with JSOC, we now limit all JSOC downloads to one connection.
  This will override all connection user settings passed to the downloader. (`#5714 <https://github.com/sunpy/sunpy/pull/5714>`__)
- Updated the ``plot`` methods on some timeseries classes to correctly label and format the time axis. (`#5720 <https://github.com/sunpy/sunpy/pull/5720>`__)
- Fixed a long-standing bug where our logger could intercept Astropy warnings in addition to SunPy warnings, and thus could conflict with Astropy's logger. (`#5722 <https://github.com/sunpy/sunpy/pull/5722>`__)
- Update asdf schemas so that references use URIs not tags as this is not
  supported by the new asdf extensions API. (`#5723 <https://github.com/sunpy/sunpy/pull/5723>`__)
- Increased the default maximum amount of records returned from HEC to 500 from 10.
  If the maximum number of records are returned, a message is shown. (`#5738 <https://github.com/sunpy/sunpy/pull/5738>`__)
- Reading a series of CDF files where at least one of them is empty no longer
  raises an error. A message for each empty file is logged at the DEBUG level. (`#5751 <https://github.com/sunpy/sunpy/pull/5751>`__)
- :func:`sunpy.map.header_helper.make_fitswcs_header` now includes a PC_ij matrix in the returned
  header if no rotation is specified. (`#5763 <https://github.com/sunpy/sunpy/pull/5763>`__)
- In the case where a map header has no PC_ij values, CROTA2 != 0, and
  CDELT1 != CDELT2, the calculation of the map rotation matrix has been fixed.
  This bug only affected maps with non-zero rotation, no PC matrix in the header,
  and un-equal scales along the two image axes. (`#5766 <https://github.com/sunpy/sunpy/pull/5766>`__)
- Maps created from :meth:`~sunpy.map.GenericMap.resample` and
  :meth:`~sunpy.map.GenericMap.superpixel` have been fixed in the case where
  the resampling was not square, and the PCi_j matrix (often a rotation matrix)
  was not a multiple of the identity matrix. When the PCi_j or CDi_j formalisms
  are used in the metadata these are now correctly modified, and the CDELT values
  are left unchanged. (`#5786 <https://github.com/sunpy/sunpy/pull/5786>`__)
- The ``__repr__`` of several ``sunpy.database`` classes have been updated to remove angular
  brackets and add equals signs. As an example, ``'<DatabaseEntry(id 3)>'`` has changed to
  ``'DatabaseEntry(id=3)'`` (`#5790 <https://github.com/sunpy/sunpy/pull/5790>`__)
- Fixed a bug when rotating a map by a matrix that is not purely a rotation.
  The likely way to inadvertently encounter this bug was when de-rotating a map with rectangular pixels that were not aligned with the coordinate axes. (`#5803 <https://github.com/sunpy/sunpy/pull/5803>`__)
- Fixed a bug where rotating a map while simultaneously scaling it could result in some of the map data being cropped out. (`#5803 <https://github.com/sunpy/sunpy/pull/5803>`__)
- Symmetric colorbar limits are no longer set on intensity images from MDI. (`#5825 <https://github.com/sunpy/sunpy/pull/5825>`__)
- Fixed plotting and peeking NORH timeseries data with ``pandas`` 1.4.0. (`#5830 <https://github.com/sunpy/sunpy/pull/5830>`__)
- In the case where ``sunpy.database.Database.fetch()`` successfully downloads only some of the search results, a ``sunpy.database.PartialFetchError`` is raised. This fixes a bug where the successful downloads would have been added to the database, but sometimes with incorrect metadata. (`#5835 <https://github.com/sunpy/sunpy/pull/5835>`__)
- When getting IRIS files from the VSO, Fido was incorrectly labelling them as XML files. (`#5868 <https://github.com/sunpy/sunpy/pull/5868>`__)
- `~sunpy.map.sources.HMIMap` now looks for ``'INSTRUME'`` instead of ``'TELESCOP'`` in order to support Helioviewer JPEG2000 versions of HMI data which do not preserve the ``'TELESCOP'`` keyword as expected in the JSOC standard. (`#5886 <https://github.com/sunpy/sunpy/pull/5886>`__)
- Fixes a bug where the ``cmap`` and ``norm`` keyword arguments were ignored when calling
  `~sunpy.map.MapSequence.plot`. (`#5889 <https://github.com/sunpy/sunpy/pull/5889>`__)
- Fix parsing of the GOES/XRS netcdf files to ignore leap seconds. (`#5915 <https://github.com/sunpy/sunpy/pull/5915>`__)
- Fixed compatibility with ``h5netcdf>0.14`` when loading GOES netcdf files. (`#5920 <https://github.com/sunpy/sunpy/pull/5920>`__)
- Fixed bugs with the rebinning and per-keV calculation for Fermi/GBM summary lightcurves (`~sunpy.timeseries.sources.GBMSummaryTimeSeries`). (`#5943 <https://github.com/sunpy/sunpy/pull/5943>`__)
- Fixed the unintentionally slow parsing of Fermi/GBM files (`~sunpy.timeseries.sources.GBMSummaryTimeSeries`). (`#5943 <https://github.com/sunpy/sunpy/pull/5943>`__)
- Fixes a bug in `~sunpy.map.sources.SJIMap` where undefined variable was
  used when parsing the wavelength.
  Also fixes the unit parsing by removing the "corrected" string from the
  ``BUNIT`` keyword as "corrected DN" cannot be parsed as a valid FITS unit. (`#5968 <https://github.com/sunpy/sunpy/pull/5968>`__)
- Fixed unit handling issue with `.GenericMap` and lowercasing the unit before it submits it to `astropy.units`. (`#5970 <https://github.com/sunpy/sunpy/pull/5970>`__)
- Fixed reading CDF files when a variable has more than 2 dimensions. If this is the case the variable will be ignored, and a user warning is provided. (`#5975 <https://github.com/sunpy/sunpy/pull/5975>`__)
- Fixed `sunpy.system_info` so it returns the extra group when an optional dependency is missing. (`#6011 <https://github.com/sunpy/sunpy/pull/6011>`__)
- Relax condition check for a HMI Synoptic map source. (`#6018 <https://github.com/sunpy/sunpy/pull/6018>`__)
- `.VSOClient` was not passing ``**kwargs`` through each download method. (`#6052 <https://github.com/sunpy/sunpy/pull/6052>`__)
- Fixed the inability to rotate images and maps with byte ordering that is different from the native byte order of the system (e.g., big-endian values on a little-endian system) for certain interpolation orders when internally using ``scikit-image``. (`#6064 <https://github.com/sunpy/sunpy/pull/6064>`__)
- Fixed a crash for dask arrays when displaying the `~sunpy.map.GenericMap` html representation. (`#6088 <https://github.com/sunpy/sunpy/pull/6088>`__)
- Constructing the color map name for a `~sunpy.map.sources.KCorMap` no longer requires the "detector" key in the metadata.
  This allows for reading files that are missing this keyword, as in the KCor JPEG2000 files. (`#6112 <https://github.com/sunpy/sunpy/pull/6112>`__)
- We now correctly pass keyword arguments in our internal FITS reader to `astropy.io.fits.open`. (`#6123 <https://github.com/sunpy/sunpy/pull/6123>`__)


Documentation
-------------

- Fixed various plotting issues with the gallery example :ref:`sphx_glr_generated_gallery_units_and_coordinates_AIA_limb_STEREO.py`. (`#5534 <https://github.com/sunpy/sunpy/pull/5534>`__)
- Improved the gallery example :ref:`sphx_glr_generated_gallery_units_and_coordinates_SDO_to_STEREO_Coordinate_Conversion.py` to better illustrate how coordinate transformations interact with submaps and coordinate plotting. (`#5534 <https://github.com/sunpy/sunpy/pull/5534>`__)
- Tidy the API Reference section of the documentation and improve the landing
  page for the docs. (`#5623 <https://github.com/sunpy/sunpy/pull/5623>`__)
- Add info about loading CDF files to the API documentation. (`#5735 <https://github.com/sunpy/sunpy/pull/5735>`__)
- Added a known issues entry about ``scikit-image`` package version pinning. (`#5865 <https://github.com/sunpy/sunpy/pull/5865>`__)
- Edited entries in the example gallery to have a consistent plotting style.
  Added said style guidelines to the example gallery page in the dev guide. (`#5870 <https://github.com/sunpy/sunpy/pull/5870>`__)
- Added the gallery example :ref:`sphx_glr_generated_gallery_map_transformations_projection_custom_origin.py`, which specifically showcases the azimuthal equidistant projection (also known as the Postel projection). (`#5961 <https://github.com/sunpy/sunpy/pull/5961>`__)
- Remove the part of the `~sunpy.map.sources.SJIMap` docstring that says
  it only works on L1 as the data work for L2 and the level checking was
  not being enforced. (`#5968 <https://github.com/sunpy/sunpy/pull/5968>`__)
- Updated the timeseries documentation to make it clear that you can pass in a numpy array. (`#6024 <https://github.com/sunpy/sunpy/pull/6024>`__)


Internal Changes
----------------

- Sped up the parsing of results from the VSO. For large queries this significantly
  reduces the time needed to perform a query to the VSO. (`#5681 <https://github.com/sunpy/sunpy/pull/5681>`__)
- `sunpy.map.GenericMap.wcs` now checks that the scale property has the correct
  units whilst constructing the WCS. (`#5682 <https://github.com/sunpy/sunpy/pull/5682>`__)
- Added `packaging <https://pypi.org/project/packaging/>`__ as a core dependency as distutils is now deprecated. (`#5713 <https://github.com/sunpy/sunpy/pull/5713>`__)
- `~sunpy.util.exceptions.SunpyWarning` is no longer a subclass of `~astropy.utils.exceptions.AstropyWarning`. (`#5722 <https://github.com/sunpy/sunpy/pull/5722>`__)
- Running the tests now requires the ``pytest-xdist`` package. By
  default tests are *not* run in parallel, but can be configured to do so
  using ``pytest-xdist`` command line options. (`#5827 <https://github.com/sunpy/sunpy/pull/5827>`__)
- Migrate the asdf infrastructure to the new style converters etc added in asdf
  2.8.0. This makes sure sunpy will be compatible with the upcoming asdf 3.0 release. (`#6057 <https://github.com/sunpy/sunpy/pull/6057>`__)
- Declare in our dependencies that we are not compatible with asdf 3.0.0 until we
  are. (`#6077 <https://github.com/sunpy/sunpy/pull/6077>`__)
- Improved performance of the code that parses dates in clients that use the
  `~sunpy.net.scraper.Scraper` to get available files. (`#6101 <https://github.com/sunpy/sunpy/pull/6101>`__)


3.1.0 (2021-10-29)
==================

Breaking Changes
----------------

- :meth:`sunpy.timeseries.sources.NOAAIndicesTimeSeries.peek` accepts ``plot_type`` as an argument instead of ``type``. (`#5200 <https://github.com/sunpy/sunpy/pull/5200>`__)
- Fill values are now set to `numpy.nan` in ``sunpy.timeseries.sources.noaa`` file
  parsers. They were previously set to a fill value of ``-1``. (`#5363 <https://github.com/sunpy/sunpy/pull/5363>`__)
- `sunpy.map.GenericMap.date` now looks for more metadata than just DATE-OBS,
  using new FITS keywords defined in version 4 of the standard.
  `sunpy.map.GenericMap.date` now returns, in order of preference:

  1. The DATE-OBS FITS keyword
  2. `~sunpy.map.GenericMap.date_average`
  3. `~sunpy.map.GenericMap.date_start`
  4. `~sunpy.map.GenericMap.date_end`
  5. The current time.

  If DATE-OBS is present alongside DATE-AVG or DATE-BEG and DATE-END, this results
  in a behaviour change to favour the new (more precisely defined) keywords.
  It is recommended
  to use `~sunpy.map.GenericMap.date_average`,
  `~sunpy.map.GenericMap.date_start`, or `~sunpy.map.GenericMap.date_end`
  instead if you need one of these specific times. (`#5449 <https://github.com/sunpy/sunpy/pull/5449>`__)
- ``sunpy.io.fits.get_header`` no longer automatically tries to add the
  WAVEUNIT keyword if it isn't present in the header. To replicate the original
  behaviour do::

    header = sunpy.io.fits.get_header(...)
    waveunit = sunpy.io.fits.extract_waveunit(header)
    if waveunit is not None:
        header['WAVEUNIT'] = waveunit

  The `sunpy.map.GenericMap.waveunit` property still uses
  ``sunpy.io.fits.extract_waveunit``` to try and get the waveunit if the
  WAVEUNIT key isn't present. (`#5501 <https://github.com/sunpy/sunpy/pull/5501>`__)
- `sunpy.map.GenericMap.wcs` no longer passes the whole ``.meta`` dictionary to
  `astropy.wcs.WCS` when constructing ``.wcs``. Instead each metadata value is
  manually taken from various map properties, which allows fixes to be made to
  the WCS without modifying the original map header. We think that
  `~sunpy.map.GenericMap.wcs` correctly sets all the keys needed for a full WCS
  header, but if you find anything missing please open an issue on the sunpy
  issue tracker. (`#5501 <https://github.com/sunpy/sunpy/pull/5501>`__)


Deprecations
------------

- ``sunpy.util.scraper.Scraper`` has been moved into `sunpy.net`, please update your imports to be ``from sunpy.net import Scraper``. (`#5364 <https://github.com/sunpy/sunpy/pull/5364>`__)
- Using "neighbour" as a resampling method in
  :func:`sunpy.image.resample.resample` is deprecated. Use "nearest" instead,
  which has the same effect. (`#5480 <https://github.com/sunpy/sunpy/pull/5480>`__)
- The `sunpy.visualization.animator` subpackage has been spun out into the
  standalone `mpl-animators <https://pypi.org/project/mpl-animators>`_ package,
  with the exception of `~sunpy.visualization.animator.MapSequenceAnimator`.
  Please update your imports to replace ``sunpy.visualization.animator`` with
  ``mpl_animators``.

  This is primarily because the ``ndcube`` package now relies on the animator
  classes as well as `sunpy`. (`#5619 <https://github.com/sunpy/sunpy/pull/5619>`__)


Removals
--------

- The deprecated ``sunpy.roi.chaincode.Chaincode`` has been removed in favour of `sunpy.net.helio.Chaincode`. (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- The deprecated ``sunpy.roi.roi`` was removed, there is no direct replacement but `astropy-regions <https://astropy-regions.readthedocs.io/en/latest/>`__ is something to consider. (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- The deprecated ``sunpy.instr`` has been removed, please use `sunkit_instruments <https://docs.sunpy.org/projects/sunkit-instruments/en/stable/>`__. (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- The deprecated ``sunpy.map.GenericMap.size`` has been removed, please use ``sunpy.map.GenericMap.data.size``. (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- The deprecated ability to read txt files from `sunpy.timeseries.sources.noaa.NOAAIndicesTimeSeries` and `sunpy.timeseries.sources.noaa.NOAAPredictIndicesTimeSeries` has been removed as the data provided by NOAA is now provided as JSON files. (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- Removed various deprecated methods on our Fido clients and responses:

  1. ``UnifiedResponse.build_table``, ``UnifiedResponse.tables``, ``UnifiedResponse.responses``, ``UnifiedResponse.get_response`` and ``UnifiedResponse.blocks`` as ``UnifiedResponse`` is now an `astropy.table.Table` that is sliceable.
  2. ``UnifiedResponse.response_block_properties`` as ``UnifiedResponse.path_format_keys`` was added as a better replacement.
  3. ``HECClient.time_query`` as you can now use ``Fido.search`` directly.
  4. ``sunpy.net.jsoc.attrs.Keys`` was not used for querying JSOC.
  5. ``sunpy.net.jsoc.JSOCClient.search_metadata`` as the functionality this provided was merged into `sunpy.net.jsoc.JSOCClient.search`.
  6. ``sunpy.net.vso.VSOClient.link`` as better search support in the client replaces this method. (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- The deprecated ``sunpy.map.GenericMap.draw_rectangle()`` has been removed, the replacement is :meth:`sunpy.map.GenericMap.draw_quadrangle` (`#5304 <https://github.com/sunpy/sunpy/pull/5304>`__)
- sunpy now errors if the unused ``.rsun`` or ``.heliographic_observer``
  attributes are set on a `~astropy.wcs.WCS`. (`#5348 <https://github.com/sunpy/sunpy/pull/5348>`__)
- Support for passing non-unit levels to :meth:`sunpy.map.GenericMap.draw_contours`
  when map data has units set has been removed, and with now raise an error. (`#5352 <https://github.com/sunpy/sunpy/pull/5352>`__)
- The ``origin`` argument to :meth:`sunpy.map.GenericMap.world_to_pixel` and
  :meth:`sunpy.map.GenericMap.pixel_to_world` has been removed. (`#5353 <https://github.com/sunpy/sunpy/pull/5353>`__)
- Support for plotting or contouring `~sunpy.map.GenericMap` on axes that are not
  `~astropy.visualization.wcsaxes.WCSAxes` has been removed. To create a
  ``WCSAxes``, use the ``projection`` argument when the axes is created, e.g.
  ``fig.add_subplot(111, projection=my_map)``. (`#5354 <https://github.com/sunpy/sunpy/pull/5354>`__)
- The following search attributes in `sunpy.net.vso.attrs` have been removed:
  ``['Time', 'Instrument', 'Wavelength', 'Source', 'Provider',
  'Level', 'Sample', 'Detector', 'Resolution', 'Physobs']``.
  Use the equivalent attribute from `sunpy.net.attrs` instead. (`#5355 <https://github.com/sunpy/sunpy/pull/5355>`__)
- The default response format from the VSO client is now a table. (`#5355 <https://github.com/sunpy/sunpy/pull/5355>`__)
- ``sunpy.net.hek.attrs.Time`` has been removed, use `sunpy.net.attrs.Time` instead. (`#5355 <https://github.com/sunpy/sunpy/pull/5355>`__)


New Features
------------

- Ensured that ``plot`` and ``peek`` will output the same figures for all `sunpy.timeseries.TimeSeries` sources. (`#5200 <https://github.com/sunpy/sunpy/pull/5200>`__)
- Added hook file and tests for using PyInstaller with sunpy. (`#5224 <https://github.com/sunpy/sunpy/pull/5224>`__)
- Allows :meth:`sunpy.map.GenericMap.draw_quadrangle` to accept pixel units as input to enable plotting boxes in the pixel space of the map, which can be different from the plot axes. (`#5275 <https://github.com/sunpy/sunpy/pull/5275>`__)
- Added the :func:`~sunpy.coordinates.propagate_with_solar_surface` context manager for transformations, which will automatically apply solar differential rotation when transforming a coordinate between frames with a change in time (``obstime``). (`#5281 <https://github.com/sunpy/sunpy/pull/5281>`__)
- Add support for parsing the observer location from a `~astropy.wcs.WCS` object
  when using the 'OBSGEO' formulation. This is the recommended way to define the
  observer location of a ground based observer. (`#5315 <https://github.com/sunpy/sunpy/pull/5315>`__)
- Added a new function, ``sunpy.visualization.draw_limb``, that draws
  the solar limb as seen from an arbitrary observer coordinate on a world
  coordinate system aware Axes. (`#5414 <https://github.com/sunpy/sunpy/pull/5414>`__)
- `sunpy.map.GenericMap.rsun_meters` now uses `sunpy.map.GenericMap.rsun_obs`
  as a fallback to calculate the assumed radius of emission if RSUN_REF metadata
  isn't present but metadata for `~sunpy.map.GenericMap.rsun_obs` is. (`#5416 <https://github.com/sunpy/sunpy/pull/5416>`__)
- Added :func:`sunpy.coordinates.utils.get_limb_coordinates` to get the solar
  limb coordinates as seen from a given observer. (`#5417 <https://github.com/sunpy/sunpy/pull/5417>`__)
- Printing the response from a `~sunpy.net.Fido` query now includes the URL where
  the data files are sourced from.

  If you develop a third-party `~sunpy.net.Fido` client, support for this can
  be automatically enabled by adding a ``info_url`` property to your
  `~sunpy.net.base_client.BaseClient` that returns a URL as a string. (`#5431 <https://github.com/sunpy/sunpy/pull/5431>`__)
- `~sunpy.timeseries.TimeSeries` can now read CDF files that conform to the
   ISTP/IACG guidelines (https://spdf.gsfc.nasa.gov/sp_use_of_cdf.html). (`#5435 <https://github.com/sunpy/sunpy/pull/5435>`__)
- The properties `~sunpy.map.GenericMap.date_start`,
  `~sunpy.map.GenericMap.date_end`, and `~sunpy.map.GenericMap.date_average` have
  been added to be drawn from the relevant FITS metadata, if present in the map
  header. (`#5449 <https://github.com/sunpy/sunpy/pull/5449>`__)
- Add default color map and normalization for `~sunpy.map.sources.HMISynopticMap`
  The default color map is 'hmimag' and the default normalization is linear between
  -1.5e-3 and +1.5e3, the expected normalization for this particular color map. (`#5464 <https://github.com/sunpy/sunpy/pull/5464>`__)
- The headers produced by :func:`~sunpy.map.header_helper.make_fitswcs_header` now include ``NAXIS``, ``NAXIS1``, and ``NAXIS2`` keywords. (`#5470 <https://github.com/sunpy/sunpy/pull/5470>`__)
- The `~astropy.wcs.WCS` instance returned by the `sunpy.map.GenericMap.wcs` property now includes the shape of the data array. (`#5470 <https://github.com/sunpy/sunpy/pull/5470>`__)
- Added the method :meth:`sunpy.map.GenericMap.reproject_to` for reprojecting a `~sunpy.map.Map` to a different WCS.
  This method requires the optional package `reproject` to be installed. (`#5470 <https://github.com/sunpy/sunpy/pull/5470>`__)
- Registered the time format ``tai_seconds`` for `astropy.time.Time` (via `~sunpy.time.TimeTaiSeconds`) to support parsing the numerical time format of TAI seconds since 1958-01-01 00:00:00.
  This format includes UTC leap seconds, and enables equivalent functionality to the ``anytim2tai`` routine in SSW. (`#5489 <https://github.com/sunpy/sunpy/pull/5489>`__)
- Added `sunpy.map.sources.WISPRMap` as a map source for WISPR on Parker Solar Probe.
  This improves the `~sunpy.map.GenericMap.name` of the map and adds correct
  information for the `~sunpy.map.GenericMap.processing_level` and
  `~sunpy.map.GenericMap.exposure_time`. (`#5502 <https://github.com/sunpy/sunpy/pull/5502>`__)
- ``sunpy.io.fits.write`` can now update the ``data`` and ``header`` of an existing HDU instance, as an alternative to creating a new instance of a specified HDU type. This adds support for writing a HDU (such as :class:`~astropy.io.fits.CompImageHDU`) initialised with non-default keyword arguments. (`#5503 <https://github.com/sunpy/sunpy/pull/5503>`__)
- Added `~sunpy.timeseries.GenericTimeSeries.observatory` to provide observatory information for the timeseries e.g. specific goes satellite number. (`#5556 <https://github.com/sunpy/sunpy/pull/5556>`__)
- :meth:`sunpy.timeseries.GenericTimeSeries.plot` and
  :meth:`sunpy.timeseries.GenericTimeSeries.peek` will now automatically label
  the y-axis if all the columns being plotted have the same units. (`#5557 <https://github.com/sunpy/sunpy/pull/5557>`__)
- :meth:`sunpy.timeseries.GenericTimeSeries.plot` and
  :meth:`sunpy.timeseries.GenericTimeSeries.peek` now have an option ``columns``
  that allows plotting a subset of the columns present. (`#5557 <https://github.com/sunpy/sunpy/pull/5557>`__)
- Added a new CDAWeb client, along with helper utilities to `sunpy.net.cdaweb`. (`#5558 <https://github.com/sunpy/sunpy/pull/5558>`__)
- Support for filtering searches with JSOC keywords has been added to ``Fido.search``. (`#5566 <https://github.com/sunpy/sunpy/pull/5566>`__)
- Added support for arithmetic operations between`~sunpy.map.GenericMap` and array-like
  objects. (`#5614 <https://github.com/sunpy/sunpy/pull/5614>`__)
- Added ``quantity`` attribute to `~sunpy.map.GenericMap` to expose the ``data``
  attribute as a `~astropy.units.Quantity` using the ``unit`` attribute. (`#5614 <https://github.com/sunpy/sunpy/pull/5614>`__)


Bug Fixes
---------

- :meth:`sunpy.map.GenericMap.superpixel` now keeps the reference coordinate of the
  WCS projection the same as the input map, and updates the reference pixel accordingly.
  This fixes inconsistencies in the input and output world coordinate systems when a
  non-linear projection is used. (`#5295 <https://github.com/sunpy/sunpy/pull/5295>`__)
- Inputs to the ``dimensions`` and ``offset`` arguments to
  :meth:`sunpy.map.GenericMap.superpixel` in units other than ``u.pix``
  (e.g. ```u.kpix``) are now handled correctly. (`#5301 <https://github.com/sunpy/sunpy/pull/5301>`__)
- Fractional inputs to the ``dimensions`` and ``offset`` arguments to
  :meth:`sunpy.map.GenericMap.superpixel` were previously rounded using `int`
  in the superpixel algorithm, but not assigned integer values in the new metadata.
  This has now been changed so the rounding is correctly reflected in the metadata. (`#5301 <https://github.com/sunpy/sunpy/pull/5301>`__)
- Remove runtime use of ``astropy.tests.helper.assert_quantity_allclose`` which
  introduces a runtime dependency on ``pytest``. (`#5305 <https://github.com/sunpy/sunpy/pull/5305>`__)
- :meth:`sunpy.map.GenericMap.resample` now keeps the reference coordinate of the
  WCS projection the same as the input map, and updates the reference pixel accordingly.
  This fixes inconsistencies in the input and output world coordinate systems when a
  non-linear projection is used. (`#5309 <https://github.com/sunpy/sunpy/pull/5309>`__)
- Fix saving `.GenericMap` to an asdf file with version 2.8.0 of the asdf package. (`#5342 <https://github.com/sunpy/sunpy/pull/5342>`__)
- When the limb is entirely visible, :meth:`sunpy.map.GenericMap.draw_limb` no
  longer plots an invisible patch for the hidden part of the limb and now returns
  `None` instead of the invisible patch. Similarly, when the limb is entirely
  invisible, no patch is drawn for the visible part and `None` is returned
  instead of the visible patch. (`#5414 <https://github.com/sunpy/sunpy/pull/5414>`__)
- :meth:`sunpy.map.GenericMap.plot` now correctly sets axis labels based on the
  coordinate system of the axes, and not the coordinate system of the map
  being plotted. This was previously only an issue if using ``autoalign=True``
  when the Map coordinate system was different to the axes coordinate system. (`#5432 <https://github.com/sunpy/sunpy/pull/5432>`__)
- :meth:`sunpy.map.GenericMap.plot` no longer adds a unit string to the axis
  labels if the axes being plotted on is a WCSAxes. For a WCSAxes, angular units
  are indicated in the tick labels, and automatically change when the zoom level
  changes from e.g. degrees to arc-minutes. This could previously lead to
  situations where the axis label units were incorrect. (`#5432 <https://github.com/sunpy/sunpy/pull/5432>`__)
- Implement automatic fallback to helioviewer mirrors if API is non-functional. (`#5440 <https://github.com/sunpy/sunpy/pull/5440>`__)
- Fixed the incorrect value for the FITS WCS ``LONPOLE`` keyword when using :func:`~sunpy.map.header_helper.make_fitswcs_header` for certain combinations of WCS projection and reference coordinate. (`#5448 <https://github.com/sunpy/sunpy/pull/5448>`__)
- The date returned by `~sunpy.map.GenericMap.date` for Solar Orbiter/EUI maps
  has been adjusted to be taken from the DATE-AVG keyword
  (the middle of the image acquisition period), instead of the DATE-OBS
  keyword (the beginning of the image acquisition period). This means the observer
  coordinate now has the correct date. (`#5462 <https://github.com/sunpy/sunpy/pull/5462>`__)
- The ``.unit`` attribute for HMI synoptic maps has been fixed. (`#5467 <https://github.com/sunpy/sunpy/pull/5467>`__)
- When "TAI" is in the date string, `sunpy.map.GenericMap.date`
  now only raises a warning if the TIMESYS keyword is present
  and different to "TAI". Previously a warning was raised all the
  time when "TAI" was in the date string. (`#5468 <https://github.com/sunpy/sunpy/pull/5468>`__)
- Fixed a bug where the property `sunpy.map.GenericMap.rsun_meters` would always internally determine the observer location, even when it is not needed, particularly for Stonyhurst heliographic maps, which have no notion of an observer.
  Thus, when working with a Stonyhurst heliographic map, a user could get an irrelevant warning message about having to assume an observer location (Earth center). (`#5478 <https://github.com/sunpy/sunpy/pull/5478>`__)
- Fixed the unintended insertion of (assumed) observer location information when accessing the property `sunpy.map.GenericMap.wcs` for Stonyhurst heliographic maps. (`#5478 <https://github.com/sunpy/sunpy/pull/5478>`__)
- Fixed an incorrect value for the FITS WCS ``LONPOLE`` keyword when using :func:`~sunpy.map.header_helper.make_fitswcs_header` for `~sunpy.coordinates.frames.Helioprojective` maps with certain values of latitude for the reference coordinate. (`#5490 <https://github.com/sunpy/sunpy/pull/5490>`__)
- A non-standard ``CROTA`` keyword included in a `sunpy.map.sources.EUIMap` FITS header is now renamed to the recommended ``CROTA2`` so a warning is no longer raised. (`#5493 <https://github.com/sunpy/sunpy/pull/5493>`__)
- The plotting x-limits of :meth:`sunpy.timeseries.sources.NOAAIndicesTimeSeries.plot`
  are now adjusted to only include finite points in the timeseries data. (`#5496 <https://github.com/sunpy/sunpy/pull/5496>`__)
- The Hinode/XRT map source now corrects the TIMESYS keyword, fixing the ``.wcs``
  property that was previously broken for Hinode/XRT maps. (`#5508 <https://github.com/sunpy/sunpy/pull/5508>`__)
- Updated `sunpy.map.CompositeMap.plot` to support the ``linestyles`` and ``colors`` arguments, in addition to the existing ``linewidths`` argument. (`#5521 <https://github.com/sunpy/sunpy/pull/5521>`__)
- Fixed a bug where rotating a `~sunpy.map.Map` could result in an extremely small shift (at the numerical-precision level) in the mapping from world coordinates to pixels. (`#5553 <https://github.com/sunpy/sunpy/pull/5553>`__)
- Fixed a bug where rotating a `~sunpy.map.Map` that is missing observation-time metadata could result in an incorrect reference coordinate. (`#5553 <https://github.com/sunpy/sunpy/pull/5553>`__)
- Fix a bug where saving a helioprojective or heliocentric coordinate to an
  asdf file didn't work due to a schema version mismatch if the observer
  location was a fully specified Stonyhurst heliographic coordinate. (`#5584 <https://github.com/sunpy/sunpy/pull/5584>`__)
- `~sunpy.map.sources.XRTMap` uppercases the ``TIMESYS`` key before checking if the
  key needs to be fixed. (`#5592 <https://github.com/sunpy/sunpy/pull/5592>`__)
- Fixed passing a URL to ``sunpy.io.read_file`` on windows. (`#5601 <https://github.com/sunpy/sunpy/pull/5601>`__)
- Fixed a bug where the ``date`` property on `~sunpy.map.sources.HMISynopticMap` returned ``None``
  if the ``DATE-OBS`` key was present. (`#5648 <https://github.com/sunpy/sunpy/pull/5648>`__)


Documentation
-------------

- Added the gallery example :ref:`sphx_glr_generated_gallery_differential_rotation_comparing_rotation_models.py` to visualize the differences between models of solar differential rotation. (`#5527 <https://github.com/sunpy/sunpy/pull/5527>`__)
- Added an example to how to save out maps as FITS files and load them back in, :ref:`sphx_glr_generated_gallery_saving_and_loading_data_genericmap_in_fits.py`. (`#5544 <https://github.com/sunpy/sunpy/pull/5544>`__)


Internal Changes
----------------

- The `~sunpy.coordinates.frames.Helioprojective` frame now has the convenience property ``angular_radius`` to return the angular radius of the Sun as seen by the observer. (`#5191 <https://github.com/sunpy/sunpy/pull/5191>`__)
- Online tests can now report back status of remote urls and will XFAIL if the remote server is unreachable. (`#5233 <https://github.com/sunpy/sunpy/pull/5233>`__)
- Re-enabled the unit test to check for coordinates consistency with JPL HORIZONS when the matching ephemeris can be specified. (`#5314 <https://github.com/sunpy/sunpy/pull/5314>`__)
- The `~sunpy.timeseries.TimeSeries` factory has been refactored to
  improve readability and maintainability of the internal code. (`#5411 <https://github.com/sunpy/sunpy/pull/5411>`__)
- `sunpy.map.GenericMap.rsun_obs` no longer emits a warning if the metadata it
  looks for is not present. Instead the standard photospheric radius is assumed
  and a log message emitted at the 'info' level. (`#5416 <https://github.com/sunpy/sunpy/pull/5416>`__)
- Nearest-neighbour and linear
  (the default for :meth:`sunpy.map.GenericMap.resample`)
  resampling have been significantly sped up. (`#5476 <https://github.com/sunpy/sunpy/pull/5476>`__)
- `sunpy.map.Map` now raises a clear error when the map is constructed if units
  of either two axes are not angular units. (`#5602 <https://github.com/sunpy/sunpy/pull/5602>`__)


3.0.1 (2021-07-03)
==================

Bug Fixes
---------

- Fixed a bug where `~sunpy.map.GenericMap` used to break with keyword arguments. (`#5392 <https://github.com/sunpy/sunpy/pull/5392>`__)
- Fixed a bug where calling :meth:`sunpy.map.GenericMap.draw_contours` on a different WCS could result in an unnecessary expansion of the plot limits. (`#5398 <https://github.com/sunpy/sunpy/pull/5398>`__)
- Fixed incorrect return values from :func:`~sunpy.map.all_corner_coords_from_map` if a rectangular map was provided. (`#5419 <https://github.com/sunpy/sunpy/pull/5419>`__)
- Do not trigger a pytest import in the asdf plugin for saving sunpy coordinate frames. (`#5429 <https://github.com/sunpy/sunpy/pull/5429>`__)
- Constructing a 2D coordinate in the `~sunpy.coordinates.frames.HeliographicCarrington` frame with ``observer='self'`` now raises an error upon creation.
  When specifying ``observer='self'``, the ``radius`` coordinate component serves as the Sun-observer distance that is necessary to fully define the Carrington heliographic coordinates. (`#5358 <https://github.com/sunpy/sunpy/pull/5358>`__)
- Fixed two bugs with handling the motion of the Sun when transforming between coordinate frames with a change in ``obstime``.
  These bugs did not affect any results if the context manager :func:`~sunpy.coordinates.transform_with_sun_center` had been used. (`#5381 <https://github.com/sunpy/sunpy/pull/5381>`__)
- Fixed a bug where the ``rsun`` frame attribute could be unintentionally reset to the default value during transformation.
  This bug primarily affected the transformation of a `~sunpy.coordinates.frames.Helioprojective` coordinate to a `~sunpy.coordinates.frames.HeliographicStonyhurst` frame. (`#5395 <https://github.com/sunpy/sunpy/pull/5395>`__)
- Fixed a bug where creating a `~sunpy.coordinates.frames.HeliographicStonyhurst` frame or a `~sunpy.coordinates.frames.HeliographicCarrington` frame from WCS information failed to make use of any specified ``rsun_ref`` value. (`#5395 <https://github.com/sunpy/sunpy/pull/5395>`__)
- `~sunpy.map.sources.SXTMap` now always returns `None` for the ``wavelength`` attribute. Previously this raised an error. (`#5401 <https://github.com/sunpy/sunpy/pull/5401>`__)


Added/Improved Documentation
----------------------------

- Simplified the "Downloading LASCO C2" gallery example by removing redundant modifications to the metadata before it is loaded by `~sunpy.map.Map`. (`#5402 <https://github.com/sunpy/sunpy/pull/5402>`__)
- Tided up the HMI synoptic map example by removing redundant code and correcting some of the comments. (`#5413 <https://github.com/sunpy/sunpy/pull/5413>`__)

3.0.0 (2021-05-14)
==================

Backwards Incompatible Changes
------------------------------

- ``sunpy.instr`` has been deprecated and will be removed in sunpy 3.1 in favour of `sunkit_instruments`.
  The code that is under ``sunpy.instr`` is imported via `sunkit_instruments` to ensure backwards comparability. (`#4526 <https://github.com/sunpy/sunpy/pull/4526>`__)
- Several `sunpy.map.GenericMap` attributes have been updated to return `None` when the relevant piece of FITS metadata is missing. These are:

  - `~sunpy.map.GenericMap.exposure_time`, previously defaulted to zero seconds.
  - `~sunpy.map.GenericMap.measurement`, previously defaulted to zero.
  - `~sunpy.map.GenericMap.waveunit`, previously defaulted to ``u.one``.
  - `~sunpy.map.GenericMap.wavelength`, previously defaulted to zero. (`#5126 <https://github.com/sunpy/sunpy/pull/5126>`__)
- `~sunpy.coordinates.frames.HeliographicStonyhurst` and `~sunpy.coordinates.frames.HeliographicCarrington` no longer automatically convert 2D input to a 3D coordinate during instantiation.
  Instead, the 2D-to-3D conversion is deferred until the coordinate is transformed to a different frame, or with a call to the method :meth:`~sunpy.coordinates.frames.BaseHeliographic.make_3d`. (`#5211 <https://github.com/sunpy/sunpy/pull/5211>`__)
- Changed URL for the `sunpy.net.dataretriever.sources.noaa.SRSClient` from "ftp://ftp.swpc.noaa.gov/pub/warehouse/" to "ftp://ftp.ngdc.noaa.gov/STP/swpc_products/daily_reports/".
  The old URL is unsupported and we expect the files will be the same but we can not say with 100% certainty. (`#5173 <https://github.com/sunpy/sunpy/pull/5173>`__)
- Changed `sunpy.net.attrs.Source` to `sunpy.net.attrs.Provider` for the `sunpy.net.dataretriever.sources.gong.GONGClient`. (`#5174 <https://github.com/sunpy/sunpy/pull/5174>`__)
- The ``rsun`` frame attribute of `~sunpy.coordinates.frames.Helioprojective` now converts any input to kilometers. (`#5211 <https://github.com/sunpy/sunpy/pull/5211>`__)
- :meth:`sunpy.map.CompositeMap.plot` now internally calls :meth:`sunpy.map.GenericMap.plot` and :meth:`sunpy.map.GenericMap.draw_contours`, which may affect the plot output of existing user code. (`#5255 <https://github.com/sunpy/sunpy/pull/5255>`__)
- Removed the ``basic_plot`` keyword argument from :meth:`sunpy.map.CompositeMap.peek` due to its unreliability. (`#5255 <https://github.com/sunpy/sunpy/pull/5255>`__)
- ``sunpy.util.sphinx.changelog`` and ``sunpy.util.towncrier`` have been removed and are now in a standalone package `sphinx-changelog <https://github.com/openastronomy/sphinx-changelog>`__. (`#5049 <https://github.com/sunpy/sunpy/pull/5049>`__)


Deprecations and Removals
-------------------------

- Deprecated ``sunpy.map.GenericMap.draw_rectangle`` in favor of :meth:`~sunpy.map.GenericMap.draw_quadrangle`. (`#5236 <https://github.com/sunpy/sunpy/pull/5236>`__)
- Using `~sunpy.map.GenericMap` plotting methods on an `~matplotlib.axes.Axes` that is not a `~astropy.visualization.wcsaxes.WCSAxes` is deprecated.
  This previously raised a warning, but is now formally deprecated, and will raise an error in sunpy 3.1. (`#5244 <https://github.com/sunpy/sunpy/pull/5244>`__)
- Deprecated ``sunpy.roi.chaincode.Chaincode`` and created a replacement at `sunpy.net.helio.Chaincode`.

  This replacement has the following changes:

  1. Added support for numpy array as an input (it was broken before).
  2. Renamed ``BoundingBox`` to ``boundingbox``
  3. Renamed ``subBoundingBox`` to ``sub_boundingbox``
  4. Now area and length raise `NotImplementedError` (`#5249 <https://github.com/sunpy/sunpy/pull/5249>`__)
- Deprecated ``sunpy.roi.roi``, as it currently has no obvious use and has never seen any real development work. (`#5249 <https://github.com/sunpy/sunpy/pull/5249>`__)


Features
--------

- :func:`sunpy.coordinates.get_horizons_coord` can now be given a start time, end time,
  and number of intervals (or interval length) to query a evenly spaced set of
  times. See the documentation string for more information and an example. (`#4698 <https://github.com/sunpy/sunpy/pull/4698>`__)
- Added :meth:`sunpy.map.GenericMap.draw_quadrangle` for drawing a quadrangle on a map.
  A quadrangle has edges that are aligned with lines of constant latitude and longitude, but these can be in a different coordinate system than that of the map. (`#4809 <https://github.com/sunpy/sunpy/pull/4809>`__)
- Added a ``longitude`` keyword argument to :func:`~sunpy.coordinates.sun.carrington_rotation_time` as an alternate way to specify a fractional Carrington rotation. (`#4879 <https://github.com/sunpy/sunpy/pull/4879>`__)
- Colorbar in `sunpy.map.GenericMap.peek` now has a unit label. (`#4930 <https://github.com/sunpy/sunpy/pull/4930>`__)
- The default axes used by ``BaseFuncAnimator.get_animation()``
  is now ``BaseFuncAnimator.axes``, instead of the currently active axes (accessed via.
  :func:`matplotlib.pyplot.gca`). The allows animations to be created on figures
  created directly using `matplotlib.figure.Figure`.

  To revert to the previous behaviour of using the current axes,
  give ``axes=plt.gca()`` to ``get_animation()``. (`#4968 <https://github.com/sunpy/sunpy/pull/4968>`__)
- Added colormaps for Solar Orbiter EUI images. These are used automatically
  when an EUI image is loaded. (`#5023 <https://github.com/sunpy/sunpy/pull/5023>`__)
- Added the ability to dynamically scale `sunpy.visualization.animator` instances.
  By specifying the ``clip_interval`` keyword, it will now clip the minimum and maximum at each slider step to the specified interval. (`#5025 <https://github.com/sunpy/sunpy/pull/5025>`__)
- Added a ``sunpy.time.timerange.TimeRange.__contains__`` method to `sunpy.time.TimeRange`
  that tests if two time ranges overlap. (`#5093 <https://github.com/sunpy/sunpy/pull/5093>`__)
- Added the ability to namespace files downloaded using `sunpy.data.data_manager.manager.DataManager` by prepending the file name with module name. (`#5111 <https://github.com/sunpy/sunpy/pull/5111>`__)
- Added a rigid rotation model to :func:`~sunpy.physics.differential_rotation.diff_rot` via ``rot_type=rigid``, where the rotation rate does not vary with latitude. (`#5132 <https://github.com/sunpy/sunpy/pull/5132>`__)
- Added a :meth:`~sunpy.map.MapSequence.save` method to `sunpy.map.MapSequence`
  that saves each map of the sequence. (`#5145 <https://github.com/sunpy/sunpy/pull/5145>`__)
- The allowable ``level`` inputs to :meth:`sunpy.map.GenericMap.contour` and
  :meth:`sunpy.map.GenericMap.draw_contours` have been consolidated. Both methods
  now accept
  - Scalars, if the map has no units
  - Quantities, if the map has units
  - Percentages (`#5154 <https://github.com/sunpy/sunpy/pull/5154>`__)
- Added support for corrected NOAA SWPC solar region summary data files. (`#5173 <https://github.com/sunpy/sunpy/pull/5173>`__)
- Updated ``sunpy.util.sysinfo.system_info`` to return all optional dependencies of sunpy. (`#5175 <https://github.com/sunpy/sunpy/pull/5175>`__)
- `sunpy.map.Map` now supports the EUI instrument on Solar Orbiter. (`#5210 <https://github.com/sunpy/sunpy/pull/5210>`__)
- `~sunpy.coordinates.frames.HeliographicStonyhurst` and `~sunpy.coordinates.frames.HeliographicCarrington` now have an ``rsun`` frame attribute to specify the radius of the Sun, which defaults to the photospheric radius defined in `sunpy.sun.constants`.
  This frame attribute is used when converting a 2D coordinate (longitude and latitude, with no specified radial distance) to a 3D coordinate by setting the radial distance to ``rsun`` (i.e., the assumption is that the coordinate is on the surface of the Sun). (`#5211 <https://github.com/sunpy/sunpy/pull/5211>`__)
- Enhanced :meth:`sunpy.map.GenericMap.draw_limb` so that the solar limb can be plotted on axes that correspond to a different map (e.g., with a different observer).
  The part of the limb that is not visible to the axes's observer because it is on the far side of the Sun is shown as dotted rather than solid. (`#5237 <https://github.com/sunpy/sunpy/pull/5237>`__)
- `~sunpy.util.MetaDict` now saves a copy of the metadata on creation, which can
  be accessed using the `~sunpy.util.MetaDict.original_meta` property. Three
  new properties have also been added to query any changes that have been made
  to metadata:

  - `~sunpy.util.MetaDict.added_items`
  - `~sunpy.util.MetaDict.removed_items`
  - `~sunpy.util.MetaDict.modified_items`

  As an example, ``my_map.meta.modified_items`` will return a dictionary mapping
  keys to their original value and current value. (`#5241 <https://github.com/sunpy/sunpy/pull/5241>`__)
- Added :func:`sunpy.map.contains_coordinate` which provides a quick way to see if a
  world coordinate is contained within the array bounds of a map. (`#5252 <https://github.com/sunpy/sunpy/pull/5252>`__)
- Added an optional keyword argument ``autoalign`` to :meth:`sunpy.map.GenericMap.plot` for plotting a map to axes that correspond to a different WCS.
  See :ref:`sphx_glr_generated_gallery_map_transformations_autoalign_aia_hmi.py`. (`#5255 <https://github.com/sunpy/sunpy/pull/5255>`__)
- :meth:`sunpy.map.CompositeMap.plot` now properly makes use of WCS information to position and orient maps when overlaying them. (`#5255 <https://github.com/sunpy/sunpy/pull/5255>`__)


Bug Fixes
---------

- Fixed the drawing methods of `sunpy.map.GenericMap` (e.g., ``~sunpy.map.GenericMap.draw_rectangle``) so that any text labels will appear in the legend. (`#5019 <https://github.com/sunpy/sunpy/pull/5019>`__)
- Fixed bug in ``sunpy.until.scraper.Scraper`` which caused URL patterns containing backslashes to be incorrectly parsed on Windows. (`#5022 <https://github.com/sunpy/sunpy/pull/5022>`__)
- Constructing a `~sunpy.util.MetaDict` is now more lenient, and accepts
  any class that inherits from `collections.abc.Mapping`. This fixes a
  regression where headers read with `astropy.io.fits` raised an error when
  passed to individual `~sunpy.map` sources. (`#5047 <https://github.com/sunpy/sunpy/pull/5047>`__)
- Added warning to :meth:`sunpy.map.GenericMap.rotate` when specified ``missing`` value is not compatible
  with the number type of the data array. (`#5051 <https://github.com/sunpy/sunpy/pull/5051>`__)
- Prevented some colormaps being accidentally modified depending on the order
  and method through which they were accessed. (`#5054 <https://github.com/sunpy/sunpy/pull/5054>`__)
- Reverted change for `sunpy.map.GenericMap.draw_limb` that made it use "add_artist" as it was changing the FOV of the plotted image. (`#5069 <https://github.com/sunpy/sunpy/pull/5069>`__)
- Fixed a bug where some `~sunpy.coordinates.metaframes.RotatedSunFrame` transformations could fail with an ``observer=None`` error. (`#5084 <https://github.com/sunpy/sunpy/pull/5084>`__)
- Fixed bug where `sunpy.data.data_manager.DataManager` would fail to recover upon deleting the sqlite database file. (`#5089 <https://github.com/sunpy/sunpy/pull/5089>`__)
- Fixed a bug where coordinate frames were considered different due to an unintended time difference during time handling at the level of numerical precision (i.e., tens of picoseconds).
  This resulted in the unexpected use of transformation machinery when transforming a coordinate to its own coordinate frame. (`#5127 <https://github.com/sunpy/sunpy/pull/5127>`__)
- Fixed a bug with failing downloads in 2010 with the `~sunpy.net.dataretriever.sources.noaa.SRSClient`. (`#5159 <https://github.com/sunpy/sunpy/pull/5159>`__)
- If the property `sunpy.map.GenericMap.rsun_obs` needs to calculate the solar angular radius from header information, it now properly uses the ``rsun_ref`` keyword if it is present and does not emit any warning. (`#5172 <https://github.com/sunpy/sunpy/pull/5172>`__)
- Added a "rsun_obs" keyword to the output of :func:`sunpy.map.header_helper.make_fitswcs_header` if the coordinate argument has a "rsun" frame attribute. (`#5177 <https://github.com/sunpy/sunpy/pull/5177>`__)
- Fixed small inaccuracies in the grid plotted by :meth:`~sunpy.map.GenericMap.draw_grid` for maps that specify a radius of the Sun that is different from the constant in `sunpy.sun.constants`. (`#5211 <https://github.com/sunpy/sunpy/pull/5211>`__)
- Fixed :meth:`sunpy.map.GenericMap.draw_contours` so that the contours from a map can be plotted on axes with a different coordinate system. (`#5239 <https://github.com/sunpy/sunpy/pull/5239>`__)
- When using the cylindrical representation of ``Heliocentric`` to work in the Heliocentric Radial coordinate frame, the ``psi`` component now goes from 0 to 360 degrees instead of -180 to 180 degrees. (`#5242 <https://github.com/sunpy/sunpy/pull/5242>`__)
- Changed ``MDIMap`` to use the "CONTENT" keyword to identify the measurement, similar to ``HMIMap``, and removed the special-case nickname. This fixes the broken title on plots. (`#5257 <https://github.com/sunpy/sunpy/pull/5257>`__)
- :func:`sunpy.coordinates.solar_frame_to_wcs_mapping` now sets the observer auxiliary
  information when a `~sunpy.coordinates.HeliographicCarrington` frame with
  ``observer='self'`` is passed. (`#5264 <https://github.com/sunpy/sunpy/pull/5264>`__)
- Calling :func:`sunpy.map.header_helper.make_fitswcs_header` with a
  `~sunpy.coordinates.HeliographicCarrington` coordinate that with ``observer='self'``
  set now correctly sets the observer information in the header. (`#5264 <https://github.com/sunpy/sunpy/pull/5264>`__)
- :meth:`sunpy.map.GenericMap.superpixel` now keeps the reference coordinate of the
  WCS projection the same as the input map, and updates the reference pixel accordingly.
  This fixes inconsistencies in the input and output world coordinate systems when a
  non-linear projection is used. (`#5295 <https://github.com/sunpy/sunpy/pull/5295>`__)
- Inputs to the ``dimensions`` and ``offset`` arguments to
  :meth:`sunpy.map.GenericMap.superpixel` in units other than ``u.pix``
  (e.g. ```u.kpix``) are now handled correctly. (`#5301 <https://github.com/sunpy/sunpy/pull/5301>`__)
- Fractional inputs to the ``dimensions`` and ``offset`` arguments to
  :meth:`sunpy.map.GenericMap.superpixel` were previously rounded using `int`
  in the superpixel algorithm, but not assigned integer values in the new metadata.
  This has now been changed so the rounding is correctly reflected in the metadata. (`#5301 <https://github.com/sunpy/sunpy/pull/5301>`__)
- Remove runtime use of ``astropy.tests.helper.assert_quantity_allclose`` which
  introduces a runtime dependency on ``pytest``. (`#5305 <https://github.com/sunpy/sunpy/pull/5305>`__)
- :meth:`sunpy.map.GenericMap.resample` now keeps the reference coordinate of the
  WCS projection the same as the input map, and updates the reference pixel accordingly.
  This fixes inconsistencies in the input and output world coordinate systems when a
  non-linear projection is used. (`#5309 <https://github.com/sunpy/sunpy/pull/5309>`__)
- Fix saving `.GenericMap` to an asdf file with version 2.8.0 of the asdf package. (`#5342 <https://github.com/sunpy/sunpy/pull/5342>`__)


Added/Improved Documentation
----------------------------

- Added a gallery example for drawing rectangles on maps. (`#4528 <https://github.com/sunpy/sunpy/pull/4528>`__)
- Added an example (:ref:`sphx_glr_generated_gallery_plotting_wcsaxes_plotting_example.py`)
  of how pixel and SkyCoords work when plotted with `~astropy.visualization.wcsaxes`. (`#4867 <https://github.com/sunpy/sunpy/pull/4867>`__)
- Added a gallery example  (:ref:`sphx_glr_generated_gallery_plotting_plotting_blank_map.py`) on how to create a blank map and mark locations. (`#5077 <https://github.com/sunpy/sunpy/pull/5077>`__)
- Added a gallery example (:ref:`sphx_glr_generated_gallery_showcase_hmi_cutout.py`)
  demonstrating how to add a HMI zoomed-in region next to a full disk HMI image. (`#5090 <https://github.com/sunpy/sunpy/pull/5090>`__)
- Updated the :ref:`sphx_glr_generated_gallery_computer_vision_techniques_mask_disk.py` example to generate the mask using :func:`sunpy.map.coordinate_is_on_solar_disk`. (`#5114 <https://github.com/sunpy/sunpy/pull/5114>`__)
- Added a gallery example (:ref:`sphx_glr_generated_gallery_map_map_segment.py`)
  demonstrating how to create a segment of a particular map from transformed coordinates. (`#5121 <https://github.com/sunpy/sunpy/pull/5121>`__)
- For the various subclasses of `~sunpy.map.GenericMap` (e.g., `~sunpy.map.sources.AIAMap`), the online documentation now shows all of the inherited attributes and methods. (`#5142 <https://github.com/sunpy/sunpy/pull/5142>`__)
- Added a documentation string to `~sunpy.map.sources.sdo.HMISynopticMap`. (`#5186 <https://github.com/sunpy/sunpy/pull/5186>`__)
- Added a new gallery example showcasing how to overlay HMI contours on an AIA image. (`#5229 <https://github.com/sunpy/sunpy/pull/5229>`__)


Trivial/Internal Changes
------------------------

- Replaced the old test runner with a new version that adds a dependency check before the test suite is run. (`#4596 <https://github.com/sunpy/sunpy/pull/4596>`__)
- The testing suite now raises a warning if the `~matplotlib.pyplot` figure stack is not empty prior to running a test, and it closes all open figures after finishing each test. (`#4969 <https://github.com/sunpy/sunpy/pull/4969>`__)
- Improved performance when moving the slider in
  ``sunpy.visualisation.animator.ArrayAnimatorWCS``. (`#4971 <https://github.com/sunpy/sunpy/pull/4971>`__)
- Added some basic logging to HEK searches, at the 'debug' logging level. (`#5020 <https://github.com/sunpy/sunpy/pull/5020>`__)
- Refactored `~sunpy.coordinates.metaframes.RotatedSunFrame` transformations for improved performance. (`#5084 <https://github.com/sunpy/sunpy/pull/5084>`__)
- Re-ordered keyword-only arguments of ``sunpy.map.GenericMap.draw_rectangle`` to match :meth:`sunpy.map.GenericMap.submap`. (`#5091 <https://github.com/sunpy/sunpy/pull/5091>`__)
- Significantly sped up calls to :func:`~sunpy.time.parse_time` for string
  arguments. This will have knock on effects, including improved performance of
  querying the VSO. (`#5108 <https://github.com/sunpy/sunpy/pull/5108>`__)
- Added tests for ``sunpy.visualization.animator.mapsequenceanimator`` and :meth:`sunpy.map.MapSequence.plot`. (`#5125 <https://github.com/sunpy/sunpy/pull/5125>`__)
- The ``CROTA`` keywords are no longer set on `sunpy.map.GenericMap.wcs`, as the
  ``PC_ij`` keywords are always set and the FITS standard says that these keywords
  must not co-exist. (`#5166 <https://github.com/sunpy/sunpy/pull/5166>`__)
- Temporarily disabled the unit test to check for coordinates consistency with JPL HORIZONS due to the inability to choose a matching ephemeris. (`#5203 <https://github.com/sunpy/sunpy/pull/5203>`__)
- :func:`~sunpy.visualization.wcsaxes_compat.wcsaxes_heliographic_overlay` now accepts ``obstime`` and ``rsun`` optional arguments.
  This function is not typically called directly by users. (`#5211 <https://github.com/sunpy/sunpy/pull/5211>`__)
- `~sunpy.map.GenericMap` plotting methods now have consistent argument
  checking for the ``axes`` argument, and will raise the same warnings
  or errors for similar ``axes`` input. (`#5223 <https://github.com/sunpy/sunpy/pull/5223>`__)
- Calling :meth:`sunpy.map.GenericMap.plot` on a
  `~astropy.visualization.wcsaxes.WCSAxes` with a different
  World Coordinate System (WCS) to the map now raises a warning,
  as the map data axes may not correctly align with the coordinate axes.
  This happens if an `~matplotlib.axes.Axes` is created with a projection
  that is a different map to the one being plotted. (`#5244 <https://github.com/sunpy/sunpy/pull/5244>`__)
- Re-enabled the unit test to check for coordinates consistency with JPL HORIZONS when the matching ephemeris can be specified. (`#5314 <https://github.com/sunpy/sunpy/pull/5314>`__)


2.1.0 (2020-02-21)
==================

Backwards Incompatible Changes
------------------------------

- Support for Python 3.6 and Numpy 1.15 has been dropped in line with
  `NEP 29 <https://numpy.org/neps/nep-0029-deprecation_policy.html>`_.
  The minimum supported version of Astropy is now 4.0, and the minimum version of scipy is now 1.2. (`#4284 <https://github.com/sunpy/sunpy/pull/4284>`__)
- Changed :func:`sunpy.coordinates.sun.B0` return type from `~astropy.coordinates.Angle`
  to `~astropy.coordinates.Latitude`. (`#4323 <https://github.com/sunpy/sunpy/pull/4323>`__)
- An error is now raised if ``vmin`` or ``vmax`` are passed to
  to `sunpy.map.GenericMap.plot` and they are already set on the map ``norm``.
  This is consistent with upcoming Matplotlib changes. (`#4328 <https://github.com/sunpy/sunpy/pull/4328>`__)
- Previously slicing the result of ``Fido.search()`` (a `~sunpy.net.fido_factory.UnifiedResponse` object) so that
  it had a length of one returned another `~sunpy.net.fido_factory.UnifiedResponse` object.
  Now it will return a `~sunpy.net.base_client.QueryResponseTable` object, which is a subclass
  of `astropy.table.Table`. (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- The ``.size`` property of a coordinate frame with no associated data will now raise
  an error instead of returning 0. (`#4577 <https://github.com/sunpy/sunpy/pull/4577>`__)
- The following `~sunpy.map.Map` methods have had support for specific positional
  arguments removed. They must now be passed as keyword arguments
  (i.e. ``m.method(keyword_arg=value)``).

  - :meth:`~sunpy.map.GenericMap.submap`: ``width``, ``height``.
  - ``sunpy.map.GenericMap.draw_rectangle``: ``width``, ``height``, ``axes``, ``top_right``.
    (`#4616 <https://github.com/sunpy/sunpy/pull/4616>`__)
- The sunpy specific attributes ``.heliographic_observer`` and ``.rsun``
  are no longer set on the `~astropy.wcs.WCS` returned by `sunpy.map.GenericMap.wcs`. (`#4620 <https://github.com/sunpy/sunpy/pull/4620>`__)
- Due to upstream changes, the parsing logic for the `~sunpy.net.helio.HECClient` now returns
  strings and not bytes for :meth:`~sunpy.net.helio.HECClient.get_table_names`. (`#4643 <https://github.com/sunpy/sunpy/pull/4643>`__)
- Reduced the selection of dependent packages installed by default via ``pip``,
  which means that some of our sub-packages will not fully import when sunpy is installed with
  ``pip install "sunpy"``.
  You can install all dependencies by specifying ``pip install "sunpy[all]"``,
  or you can install sub-package-specific dependencies by specifying, e.g.,
  ``[map]`` or ``[timeseries]``. (`#4662 <https://github.com/sunpy/sunpy/pull/4662>`__)
- The class inheritance for `~sunpy.coordinates.metaframes.RotatedSunFrame` and the frames it
  creates has been changed in order to stop depending on unsupported behavior in the underlying machinery.
  The return values for some :func:`isinstance`/:func:`issubclass` calls will be different,
  but the API for `~sunpy.coordinates.metaframes.RotatedSunFrame` is otherwise unchanged. (`#4691 <https://github.com/sunpy/sunpy/pull/4691>`__)
- Fix a bug in `~sunpy.map.GenericMap.submap` where only the top right and bottom
  left coordinates of the input rectangle in world coordinates were considered
  when calculating the pixel bounding box. All four corners are once again taken
  into account now, meaning that `~sunpy.map.GenericMap.submap` correctly returns
  the smallest pixel box which contains all four corners of the input rectangle.

  To revert to the previous 2.0.0 behaviour, first convert the top right and bottom
  left coordinates to pixel space before calling submap with::

      top_right = smap.wcs.world_to_pixel(top_right) * u.pix
      bottom_left = smap.wcs.world_to_pixel(bottom_left) * u.pix
      smap.submap(bottom_left=bottom_left, top_right=top_right)

  This will define the rectangle in pixel space. (`#4727 <https://github.com/sunpy/sunpy/pull/4727>`__)
- VSO results where the size was ``-1`` (missing data) now return ``None`` rather
  than ``-1`` to be consistent with other missing data in the VSO results. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- All result objects contained within the results of a ``Fido.search()`` (a
  `~sunpy.net.fido_factory.UnifiedResponse` object) are now
  `~sunpy.net.base_client.QueryResponseTable` objects (or subclasses thereof).
  These objects are subclasses of `astropy.table.Table` and can therefore be
  filtered and inspected as tabular objects, and the modified tables can be passed
  to ``Fido.fetch``.

  This, while a breaking change for anyone accessing these response objects
  directly, will hopefully make working with ``Fido`` search results much easier. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- Results from the `~sunpy.net.dataretriever.NOAAIndicesClient` and the
  `~sunpy.net.dataretriever.NOAAPredictClient` no longer has ``Start Time`` or
  ``End Time`` in their results table as the results returned from the client are
  not dependent upon the time parameter of a search. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- The ``sunpy.net.vso.QueryResponse.search`` method has been removed as it has not
  worked since the 1.0 release of sunpy. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- The ``sunpy.net.hek.hek.HEKColumn`` class has been removed, the ``HEKTable`` class
  now uses the standard `astropy.table.Column` class. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- The keys used to format file paths in ``Fido.fetch`` have changed. They are now
  more standardised across all the clients, as they are all extracted from the
  names of the columns in the results table.

  For results from the VSO the keys are no longer separated with ``.``, and are
  based on the displayed column names. For results from the ``dataretriever``
  clients the only main change is that the keys are now lower case, where they
  were capitalized before. You can use the ``.sunpy.net.fido_factory.UnifiedResponse.path_format_keys``
  method to see all the possible keys for a particular search. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- The time returned from :func:`~sunpy.coordinates.sun.carrington_rotation_number`
  has been changed from the TT scale to the more common UTC scale. To undo this change,
  use ``time_out = time_out.tt`` on the outputted time. (`#4819 <https://github.com/sunpy/sunpy/pull/4819>`__)
- `~.BaseQueryResponse.response_block_properties` has been renamed to
  ``.BaseQueryResponse.path_format_keys``, on the return objects from all
  ``search()`` methods on all clients and from ``Fido.search()``. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)

Removals
--------
- Removed deprecated functions:

  - ``sunpy.coordinates.frames.Helioprojective.calculate_distance``, alternative
    is `sunpy.coordinates.frames.Helioprojective.make_3d`.
  - ``sunpy.image.coalignment.repair_image_nonfinite`` - if you wish to repair the image,
    this has to be done manually before calling the various ``sunpy.image.coalignment`` functions.
  - The ``repair_nonfinite`` keyword argument to ``calculate_shift`` and  ``calculate_match_template_shift``
    has been removed.
  - ``sunpy.instr.lyra.download_lytaf_database`` - this just downloaded the file
    at ``http://proba2.oma.be/lyra/data/lytaf/annotation_ppt.db``, which can be done manually.
  - ``sunpy.util.net.check_download_file``, no alternative.
  - ``sunpy.visualization.animator.ImageAnimatorWCS``, alternative is
    ``sunpy.visualization.animator.ArrayAnimatorWCS``. (`#4350 <https://github.com/sunpy/sunpy/pull/4350>`__)

- Removed deprecated function ``sunpy.instr.aia.aiaprep``.
  Alternative is `~aiapy.calibrate.register` for converting AIA
  images from level 1 to level 1.5. (`#4485 <https://github.com/sunpy/sunpy/pull/4485>`__)
- ``sunpy.cm`` has been removed. All of the functionality in this module can
  now be found in `sunpy.visualization.colormaps`. (`#4488 <https://github.com/sunpy/sunpy/pull/4488>`__)
- ``sunpy.test.hash`` has been removed, the functionality has been moved into the
  `pytest-mpl <https://github.com/matplotlib/pytest-mpl>`__ package. (`#4605 <https://github.com/sunpy/sunpy/pull/4605>`__)
- ``sunpy.util.multimethod`` has been removed. (`#4614 <https://github.com/sunpy/sunpy/pull/4614>`__)
- The ``lytaf_path`` argument (which previously did nothing) has been removed from
  - ``sunpy.instr.lyra.remove_lytaf_events_from_timeseries``
  - ``sunpy.instr.lyra.get_lytaf_events``
  - ``sunpy.instr.lyra.get_lytaf_event_types`` (`#4615 <https://github.com/sunpy/sunpy/pull/4615>`__)

Deprecations
------------
- Deprecated ``sunpy.net.vso.attrs.Source`` and ``sunpy.net.vso.attrs.Provider``.
  They are now `sunpy.net.attrs.Source` and `sunpy.net.attrs.Provider` respectively.
  (`#4321 <https://github.com/sunpy/sunpy/pull/4321>`__)
- Deprecated the use of the ``sunpy.map.GenericMap.size`` property,
  use ``sunpy.map.Map.data.size`` instead. (`#4338 <https://github.com/sunpy/sunpy/pull/4338>`__)
- ``sunpy.net.helio.HECClient.time_query`` is deprecated, `~sunpy.net.helio.HECClient.search`
  is the replacement. (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- ``sunpy.net.jsoc.attrs.Keys`` is deprecated; all fields are returned by default and can be filtered post search. (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- ``sunpy.net.hek.attrs.Time`` is deprecated; `~sunpy.net.attrs.Time` should be used instead. (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- Support for :func:`sunpy.coordinates.wcs_utils.solar_wcs_frame_mapping` to
  use the ``.heliographic_observer`` and ``.rsun`` attributes on a
  `~astropy.wcs.WCS` is deprecated. (`#4620 <https://github.com/sunpy/sunpy/pull/4620>`__)
- The ``origin`` argument to `sunpy.map.GenericMap.pixel_to_world` and
  `sunpy.map.GenericMap.world_to_pixel` is deprecated.

  - If passing ``0``, not using the ``origin`` argument will have the same effect.
  - If passing ``1``, manually subtract 1 pixel from the input to ``pixel_to_world``,
    or manually add 1 pixel to the output of ``world_to_pixel``, and do not use the
    ``origin`` argument. (`#4700 <https://github.com/sunpy/sunpy/pull/4700>`__)
- The ``.VSOClient.link`` method is deprecated as it is no longer used. (`#4789 <https://github.com/sunpy/sunpy/pull/4789>`__)
- The ``.UnifiedResponse.get_response``, ``.UnifiedResponse.tables`` and
  ``.UnifiedResponse.responses`` attributes of ``.UnifiedResponse`` have been
  deprecated as they are no longer needed now the object returns the table
  objects it contains when sliced. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- :meth:`sunpy.net.vso.VSOClient.search` has a new keyword argument
  ``response_type=`` which controls the return type from the ``search()`` method.
  In sunpy 2.1 and 3.0 it will default to the ``"legacy"`` response format, in
  3.1 it will default to the new ``"table"`` response format, and the
  ``"legacy"`` format may be deprecated and removed at a later date.

  Searches made with ``Fido`` will use the new ``"table"`` response format, so
  this only affects users interacting with the ``VSOClient`` object directly. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)


Features
--------


- For :meth:`sunpy.map.GenericMap.quicklook` and :meth:`sunpy.map.MapSequence.quicklook` (also used for the HTML representation shown in Jupyter notebooks), the histogram is now shaded corresponding to the colormap of the plotted image.
  Clicking on the histogram will toggle an alternate version of the histogram. (`#4931 <https://github.com/sunpy/sunpy/pull/4931>`__)
- Add an ``SRS_TABLE`` file to the sample data, and use it in the magnetogram
  plotting example. (`#4993 <https://github.com/sunpy/sunpy/pull/4993>`__)
- Added a `sunpy.map.GenericMap.contour()` method to find the contours on a map. (`#3909 <https://github.com/sunpy/sunpy/pull/3909>`__)
- Added a context manager (:meth:`~sunpy.coordinates.frames.Helioprojective.assume_spherical_screen`)
  to interpret `~sunpy.coordinates.frames.Helioprojective` coordinates as being on
  the inside of a spherical screen instead of on the surface of the Sun. (`#4003 <https://github.com/sunpy/sunpy/pull/4003>`__)
- Added `sunpy.map.sources.HMISynopticMap` for handling the Synoptic maps from HMI. (`#4053 <https://github.com/sunpy/sunpy/pull/4053>`__)
- Added a `~sunpy.map.sources.MDISynopticMap` map source class. (`#4054 <https://github.com/sunpy/sunpy/pull/4054>`__)
- Created `~sunpy.net.dataretriever.GONGClient` for accessing magnetogram synoptic map archives of NSO-GONG. (`#4055 <https://github.com/sunpy/sunpy/pull/4055>`__)
- All coordinate frames will now show the velocity if it exists in the underlying data. (`#4102 <https://github.com/sunpy/sunpy/pull/4102>`__)
- The ephemeris functions :func:`~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst()`, :func:`~sunpy.coordinates.ephemeris.get_earth()`, and :func:`~sunpy.coordinates.ephemeris.get_horizons_coord()` can now optionally return the body's velocity as part of the output coordinate. (`#4102 <https://github.com/sunpy/sunpy/pull/4102>`__)
- `~sunpy.util.metadata.MetaDict` now maintains coherence between its keys and their corresponding keycomments. Calling ``del`` on a ``MetaDict`` object key is now case-insensitive. (`#4129 <https://github.com/sunpy/sunpy/pull/4129>`__)
- Allow ``sunpy.visualization.animator.ArrayAnimatorWCS`` to disable ticks for
  a coordinate, by setting ``ticks: False`` in the ``coord_params`` dictionary. (`#4270 <https://github.com/sunpy/sunpy/pull/4270>`__)
- Added a ``show()`` method for `~sunpy.net.base_client.BaseQueryResponse` which returns `~astropy.table.Table` with specified columns for the Query Response. (`#4309 <https://github.com/sunpy/sunpy/pull/4309>`__)
- Added ``_extract_files_meta`` method in ``sunpy.util.scraper.Scraper`` which allows scraper to extract metadata from the file URLs retrieved for a given time range. (`#4313 <https://github.com/sunpy/sunpy/pull/4313>`__)
- Refactoring of `~sunpy.net.dataretriever` which adds these capabilities to `~sunpy.net.dataretriever.QueryResponse`:

  - Any ``attr`` shall not be defaulted to a hard-coded value in all subclasses of `~sunpy.net.dataretriever.GenericClient`; thus records for all possible ``attrs`` shall be returned if it is not specified in the query.
  - `~sunpy.net.dataretriever.QueryResponse` can now show more columns; thus all metadata extractable from matching file URLs shall be shown and for a client, non-supported ``attrs`` shall not be shown in the response tables. (`#4321 <https://github.com/sunpy/sunpy/pull/4321>`__)
- New class attributes added to `~sunpy.net.dataretriever.GenericClient`:

  - ``baseurl`` and ``pattern`` which are required to define a new simple client.
  - ``optional`` and ``required`` which are a ``set`` of optional and required `~sunpy.net.attrs` respectively; which generalizes :meth:`~sunpy.net.dataretriever.GenericClient._can_handle_query`. (`#4321 <https://github.com/sunpy/sunpy/pull/4321>`__)
- Additions in ``sunpy.util.scraper`` to support the refactoring of `~sunpy.net.dataretriever.GenericClient`:
  - ``sunpy.util.scraper.Scraper.findDatewith_extractor`` that parses the url using extractor to return its start time.
  - A ``matcher`` in ``sunpy.util.scraper.Scraper._extract_files_meta`` which validates the extracted metadata by using the dictionary returned from :meth:`~sunpy.net.dataretriever.GenericClient._get_match_dict`. (`#4321 <https://github.com/sunpy/sunpy/pull/4321>`__)
- Added methods :meth:`~sunpy.net.dataretriever.GenericClient.pre_search_hook` and :meth:`~sunpy.net.dataretriever.GenericClient.post_search_hook` which helps to translate the attrs for scraper before and after the search respectively. (`#4321 <https://github.com/sunpy/sunpy/pull/4321>`__)
- ``sunpy.timeseries.sources.RHESSISummaryTimeSeries.peek`` has had the following minor
  changes:

  - Colors from the default matplotlib color cycle are now used (but the colors remain qualitatively the same)
  - The default matplotlib linewidth is now used
  - It is now possible to pass in a user specified linewidth
  - Seconds have been added to the x-axis labels (previously it was just hours and minutes) (`#4326 <https://github.com/sunpy/sunpy/pull/4326>`__)
- `~sunpy.net.helio.hec.HECClient` and  `~sunpy.net.hek.hek.HEKClient` now inherit `~sunpy.net.base_client.BaseClient` which makes them compatible with the `~sunpy.net.fido_factory.UnifiedDownloaderFactory` (``Fido``). (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- `~sunpy.net.helio.attrs.MaxRecords` and `~sunpy.net.helio.attrs.TableName` added as "attrs" for HELIO searches. (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- Add the ability to download new GOES 16 & 17 data alongside the reprocessed GOES 13, 14 and 15 data via the GOES-XRS Fido client. (`#4394 <https://github.com/sunpy/sunpy/pull/4394>`__)
- `sunpy.net.jsoc.JSOCClient.request_data` now support additional parameter "method" which allows user to download staged data as single .tar file. (`#4405 <https://github.com/sunpy/sunpy/pull/4405>`__)
- Added ``sunpy.util.get_timerange_from_exdict`` which finds time range for a URL using its metadata.
  Added ``sunpy.util.scraper.Scraper.isvalid_time`` that checks whether the file corresponds to a desired time range. (`#4419 <https://github.com/sunpy/sunpy/pull/4419>`__)
- Colormap data has been moved to individual .csv files in the
  :file:`sunpy/visualization/colormaps/data` directory. (`#4433 <https://github.com/sunpy/sunpy/pull/4433>`__)
- Added `~sunpy.coordinates.utils.solar_angle_equivalency` to convert between a physical distance on the Sun (e.g., km) to an angular separation as seen by an observer (e.g., arcsec). (`#4443 <https://github.com/sunpy/sunpy/pull/4443>`__)
- `sunpy.map.Map` instances now have their ``.unit`` attribute set from the
  ``'BUNIT'`` FITS keyword. If the keyword cannot be parsed, or is not present
  the unit is set to `None`. (`#4451 <https://github.com/sunpy/sunpy/pull/4451>`__)
- The `sunpy.map.GenericMap.wcs` property is now cached, and will be recomputed
  only if changes are made to the map metadata. This improves performance of a
  number of places in the code base, and only one warning will now be raised
  about WCS fixes for a given set of metadata (as opposed to a warning each time
  ``.wcs`` is accessed) (`#4467 <https://github.com/sunpy/sunpy/pull/4467>`__)
- Extended :meth:`~sunpy.timeseries.GenericTimeSeries.concatenate` and
  :meth:`~sunpy.timeseries.TimeSeriesMetaData.concatenate` to allow iterables. (`#4499 <https://github.com/sunpy/sunpy/pull/4499>`__)
- Enable `~sunpy.coordinates.metaframes.RotatedSunFrame` to work with non-SunPy frames (e.g., `~astropy.coordinates.HeliocentricMeanEcliptic`). (`#4577 <https://github.com/sunpy/sunpy/pull/4577>`__)
- Add support for `pathlib.Path` objects to be passed to `sunpy.timeseries.TimeSeries`. (`#4589 <https://github.com/sunpy/sunpy/pull/4589>`__)
- Add support for GOES XRS netcdf files to be read as a `sunpy.timeseries.sources.XRSTimeSeries`. (`#4592 <https://github.com/sunpy/sunpy/pull/4592>`__)
- Add `~sunpy.net.jsoc.attrs.Cutout` attr for requesting cutouts
  from JSOC via `~sunpy.net.jsoc.JSOCClient` and ``Fido``. (`#4595 <https://github.com/sunpy/sunpy/pull/4595>`__)
- sunpy now sets auxiliary parameters on `sunpy.map.GenericMap.wcs` using the
  `astropy.wcs.Wcsprm.aux` attribute. This stores observer information, along with
  the reference solar radius if present. (`#4620 <https://github.com/sunpy/sunpy/pull/4620>`__)
- The `~sunpy.coordinates.frames.HeliographicCarrington` frame now accepts the specification of ``observer='self'`` to indicate that the coordinate itself is also the observer for the coordinate frame.
  This functionality greatly simplifies working with locations of observatories that are provided in Carrington coordinates. (`#4659 <https://github.com/sunpy/sunpy/pull/4659>`__)
- Add two new colormaps (``rhessi`` and ``std_gamma_2``) that are used for plotting RHESSI maps. (`#4665 <https://github.com/sunpy/sunpy/pull/4665>`__)
- If either 'CTYPE1' or 'CTYPE2' are not present in map metadata, sunpy now assumes
  they are 'HPLN-TAN' and 'HPLT-TAN' (previously it assumed 'HPLN-   ' and 'HPLT-   ').
  In addition, a warning is also now raised when this assumption is made. (`#4702 <https://github.com/sunpy/sunpy/pull/4702>`__)
- Added a new `~sunpy.map.all_corner_coords_from_map` function to get the
  coordinates of all the pixel corners in a `~sunpy.map.GenericMap`. (`#4776 <https://github.com/sunpy/sunpy/pull/4776>`__)
- Added support for "%Y/%m/%dT%H:%M" to :func:`sunpy.time.parse_time`. (`#4791 <https://github.com/sunpy/sunpy/pull/4791>`__)
- Added the STEREO EUVI instrument specific colormaps called" 'euvi171', 'euvi195', 'euvi284', 'euvi304'. (`#4822 <https://github.com/sunpy/sunpy/pull/4822>`__)


Bug Fixes
---------

- `sunpy.map.GenericMap.date` now has its time scale set from the 'TIMESYS' FITS keyword,
  if it is present. If it isn't present the time scale defaults to 'UTC', which is unchanged
  default behaviour, so this change will only affect maps with a 'TIMESYS' keyword
  that is not set to 'UTC'. (`#4881 <https://github.com/sunpy/sunpy/pull/4881>`__)
- Fixed the `sunpy.net.dataretriever.sources.noaa.SRSClient` which silently failed to download the SRS files when the tarball for the previous years did not exist.
  Client now actually searches for the tarballs and srs files on the ftp archive before returning them as results. (`#4904 <https://github.com/sunpy/sunpy/pull/4904>`__)
- No longer is the WAVEUNIT keyword injected into a data source if it is missing from the file's metadata. (`#4926 <https://github.com/sunpy/sunpy/pull/4926>`__)
- Map sources no longer overwrite FITS metadata keywords if they are present in
  the original metadata. The particular map sources that have been fixed are
  `~sunpy.map.sources.SJIMap`, `~sunpy.map.sources.KCorMap`, `~sunpy.map.sources.RHESSIMap`,
  `~sunpy.map.sources.EITMap`, `~sunpy.map.sources.EUVIMap`, `~sunpy.map.sources.SXTMap`. (`#4926 <https://github.com/sunpy/sunpy/pull/4926>`__)
- Fixed a handling bug in ``sunpy.map.GenericMap.draw_rectangle`` when the rectangle is specified in a different coordinate frame than that of the map.
  A couple of other minor bugs in ``sunpy.map.GenericMap.draw_rectangle`` were also fixed. (`#4929 <https://github.com/sunpy/sunpy/pull/4929>`__)
- Improved error message from ``sunpy.net.Fido.fetch()`` when no email has been supplied for JSOC data. (`#4950 <https://github.com/sunpy/sunpy/pull/4950>`__)
- Fixed a bug when transforming from `~sunpy.coordinates.metaframes.RotatedSunFrame` to another frame at a different observation time that resulted in small inaccuracies.
  The translational motion of the Sun was not being handled correctly. (`#4979 <https://github.com/sunpy/sunpy/pull/4979>`__)
- Fixed two bugs with :func:`~sunpy.physics.differential_rotation.differential_rotate` and :func:`~sunpy.physics.differential_rotation.solar_rotate_coordinate` that resulted in significant inaccuracies.
  Both functions now ignore the translational motion of the Sun. (`#4979 <https://github.com/sunpy/sunpy/pull/4979>`__)
- The ability to to filter search results from the `~sunpy.net.vso.VSOClient` was broken.
  This has now been restored. (`#4011 <https://github.com/sunpy/sunpy/pull/4011>`__)
- Fixed a bug where transformation errors were not getting raised in some situations when a coordinate frame had ``obstime`` set to the default value of ``None`` and `~astropy.coordinates.SkyCoord` was not being used.
  Users are recommended to use `~astropy.coordinates.SkyCoord` to manage coordinate transformations unless they have a specific reason not to. (`#4267 <https://github.com/sunpy/sunpy/pull/4267>`__)
- Fixed a bug in `~sunpy.net.dataretriever.sources.goes.XRSClient._get_url_for_timerange` which returned incorrect URLs
  because of not using ``**kwargs`` in the client's ``_get_overlap_urls()`` method. (`#4288 <https://github.com/sunpy/sunpy/pull/4288>`__)
- Data products from `~sunpy.net.dataretriever.NOAAIndicesClient` and
  `~sunpy.net.dataretriever.NOAAPredictClient` have been updated to download
  new JSON files. The old text files which the data used to come in no longer
  exist. The new JSON files for `~sunpy.net.dataretriever.NOAAIndicesClient`
  now do not have the following columns:
  - Geomagnetic Observed and Smoothed
  - Sunspot Numbers Ratio (RI/SW)

  Both `sunpy.timeseries.sources.noaa.NOAAIndicesTimeSeries` and
  `sunpy.timeseries.sources.noaa.NOAAPredictIndicesTimeSeries` have been updated to
  support the new JSON files. Loading the old text files is still supported,
  but support for this will be removed in a future version of sunpy. (`#4340 <https://github.com/sunpy/sunpy/pull/4340>`__)
- Fixed a bug due to which ``sunpy.net.helio.parser.wsdl_retriever`` ignored previously discovered Taverna links. (`#4358 <https://github.com/sunpy/sunpy/pull/4358>`__)
- The flare class labels in GOES ``peek()`` plots are now drawn at the center of
  the flare classes. Previously they were (ambiguously) drawn on the boundaries. (`#4364 <https://github.com/sunpy/sunpy/pull/4364>`__)
- `sunpy.map.GenericMap.rsun_obs` no longer assumes the observer is at Earth if
   ``rsun_obs`` was not present in the map metadata. The sun-observer
   distance is now taken directly from the observer coordinate. If the observer
   coordinate is not present, this defaults to the Earth, retaining previous
   behaviour. (`#4375 <https://github.com/sunpy/sunpy/pull/4375>`__)
- Nanosecond precision is now retained when using `~sunpy.time.parse_time` with
  a `~pandas.Timestamp`. (`#4409 <https://github.com/sunpy/sunpy/pull/4409>`__)
- Fixed a bug where SunPy could not be successfully imported if the default text encoding of the running environment was unable to handle non-ASCII characters. (`#4422 <https://github.com/sunpy/sunpy/pull/4422>`__)
- `sunpy.net.dataretriever.sources.noaa.SRSClient` now correctly returns zero
  results for queries in the future or before 1996, which is when data is first
  available. (`#4432 <https://github.com/sunpy/sunpy/pull/4432>`__)
- Fixes issue where NAXISn is not updated after invoking :meth:`.GenericMap.resample` (`#4445 <https://github.com/sunpy/sunpy/pull/4445>`__)
- The floating point precision of input to `sunpy.image.transform.affine_transform`
  is now preserved. Previously all input was cast to `numpy.float64`, which could
  cause large increases in memory use for 32 bit data. (`#4452 <https://github.com/sunpy/sunpy/pull/4452>`__)
- Fixed :func:`~sunpy.image.transform.affine_transform` to scale images to [0, 1] before
  passing them to :func:`skimage.transform.warp` and later rescale them back. (`#4477 <https://github.com/sunpy/sunpy/pull/4477>`__)
- Several ``warnings.simplefilter('always', Warning)`` warning filters in
  `sunpy.timeseries` have been removed. (`#4511 <https://github.com/sunpy/sunpy/pull/4511>`__)
- All calculations of the angular radius of the Sun now use the same underlying code with the accurate calculation.
  The previous inaccuracy was a relative error of ~0.001% (0.01 arcseconds) for an observer at 1 AU, but could be as large as ~0.5% for Parker Solar Probe perihelia. (`#4524 <https://github.com/sunpy/sunpy/pull/4524>`__)
- Fixed an issue in :meth:`sunpy.time.TimeRange.get_dates` where the function would return the wrong number of days if less than 24 hours had passed (`#4529 <https://github.com/sunpy/sunpy/pull/4529>`__)
- Several functions in `sunpy.map` now properly check if the provided coordinate is in the expected `~sunpy.coordinates.frames.Helioprojective` frame. (`#4552 <https://github.com/sunpy/sunpy/pull/4552>`__)
- Fixes a bug which occurs in setting the ``ylims`` by ``sunpy.visualization.animator.line.LineAnimator`` when there are non-finite values in the data array to be animated. (`#4554 <https://github.com/sunpy/sunpy/pull/4554>`__)
- Clear rotation metadata for SOHO/LASCO Helioviewer JPEG2000 images, as they are already rotated correctly. (`#4561 <https://github.com/sunpy/sunpy/pull/4561>`__)
- The ``max_conn`` argument to ``Fido.fetch()`` is now correctly respected by
  the JSOC client. Previously the JSOC client would default to 4 connections no
  matter what the value passed to ``Fido.fetch()`` was. (`#4567 <https://github.com/sunpy/sunpy/pull/4567>`__)
- :func:`sunpy.time.parse_time` now correctly parses lists of time strings that
  have one of the built in sunpy time formats. (`#4590 <https://github.com/sunpy/sunpy/pull/4590>`__)
- Fixes the SRSClient to search for files of correct queried time and now allows a path keyword to be downloaded in fetch. (`#4600 <https://github.com/sunpy/sunpy/pull/4600>`__)
- Fixed ``sunpy.net.helio.parser.wsdl_retriever``, which previously
  ignored discovered Taverna links. (`#4601 <https://github.com/sunpy/sunpy/pull/4601>`__)
- The transformations between `~astropy.coordinates.HCRS` and `~sunpy.coordinates.frames.HeliographicStonyhurst` have been re-implemented to enable the proper transformations of velocities.
  All ephemeris functions (e.g., :func:`~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst`) now return properly calculated velocities when ``include_velocity=True`` is specified. (`#4613 <https://github.com/sunpy/sunpy/pull/4613>`__)
- The maximum number of connections opened by the JSOC downloader has been reduced
  from 4 to 2. This should prevent downloads of large numbers of files crashing. (`#4624 <https://github.com/sunpy/sunpy/pull/4624>`__)
- Fixed a significant performance bug that affected all coordinate transformations.
  Transformations have been sped up by a factor a few. (`#4663 <https://github.com/sunpy/sunpy/pull/4663>`__)
- Fixed a bug with the mapping of a WCS header to a coordinate frame if the observer location is provided in Carrington coordinates. (`#4669 <https://github.com/sunpy/sunpy/pull/4669>`__)
- ``sunpy.io.fits.header_to_fits`` now excludes any keys that have associated NaN
  values, as these are not valid in a FITS header, and throws a warning if this
  happens. (`#4676 <https://github.com/sunpy/sunpy/pull/4676>`__)
- Fixed an assumption in `sunpy.map.GenericMap.pixel_to_world` that the first
  data axis is longitude, and the second is latitude. This will affect you if
  you are using data where the x/y axes are latitude/longitude, and now returns
  correct values in methods and properties that call ``pixel_to_world``,
  such as ``bottom_left_coord``, ``top_right_coord``, ``center``. (`#4700 <https://github.com/sunpy/sunpy/pull/4700>`__)
- Added a warning when a 2D `~sunpy.coordinates.frames.Helioprojective` coordinate is upgraded to a 3D coordinate and the number type is lower precision than the native Python float.
  This 2D->3D upgrade is performed internally when transforming a 2D `~sunpy.coordinates.frames.Helioprojective` coordinate to any other coordinate frame. (`#4724 <https://github.com/sunpy/sunpy/pull/4724>`__)
- All columns from a :meth:`sunpy.net.vso.vso.VSOClient.search` will now be shown. (`#4788 <https://github.com/sunpy/sunpy/pull/4788>`__)
- The search results object returned from ``Fido.search``
  (`~sunpy.net.fido_factory.UnifiedResponse`) now correctly counts all results in
  it's `~sunpy.net.fido_factory.UnifiedResponse.file_num` property. Note that
  because some ``Fido`` clients now return metadata only results, this is really
  the number of records and does not always correspond to the number of files
  that would be downloaded. (`#4798 <https://github.com/sunpy/sunpy/pull/4798>`__)
- Improved the file processing logic for EVE L0CS files, which may have fixed a
  bug where the first line of data was parsed incorrectly. (`#4805 <https://github.com/sunpy/sunpy/pull/4805>`__)
- Fixing the ``CROTA`` meta keyword in EUVI FITS to ``CROTAn`` standard. (`#4846 <https://github.com/sunpy/sunpy/pull/4846>`__)


Added/Improved Documentation
----------------------------

- Added a developer guide for writing a new ``Fido`` client. (`#4387 <https://github.com/sunpy/sunpy/pull/4387>`__)
- Added an example of how to use Matplotlib's axes range functionality when plotting a Map with WCSAxes. (`#4792 <https://github.com/sunpy/sunpy/pull/4792>`__)
- Add links to Thompson 2006 paper on solar coordinates to synoptic map example. (`#3549 <https://github.com/sunpy/sunpy/pull/3549>`__)
- Clarified the meaning of ``.bottom_left_coord`` and ``.top_right_coord`` in
  `sunpy.map.GenericMap`. (`#3706 <https://github.com/sunpy/sunpy/pull/3706>`__)
- Added a list of possible signatures to
  `sunpy.timeseries.metadata.TimeSeriesMetaData`. (`#3709 <https://github.com/sunpy/sunpy/pull/3709>`__)
- Added `sunpy.data.manager`, `sunpy.data.cache`, `sunpy.net.Fido`, `sunpy.map.Map`,
  and `sunpy.timeseries.TimeSeries` to the docs. (`#4098 <https://github.com/sunpy/sunpy/pull/4098>`__)
- Clarified spline option for `sunpy.map.GenericMap.resample`. (`#4136 <https://github.com/sunpy/sunpy/pull/4136>`__)
- Updated the gallery example :ref:`sphx_glr_generated_gallery_plotting_solar_cycle_example.py` to retrieve data using `~sunpy.net.Fido`. (`#4169 <https://github.com/sunpy/sunpy/pull/4169>`__)
- Fixed example usage of ``sunpy.io.fits.read`` to account for the fact that it returns a list
  of data-header pairs rather than the data-header pairs directly. (`#4183 <https://github.com/sunpy/sunpy/pull/4183>`__)
- Added example of how to create a `sunpy.map.GenericMap` from observations in RA-DEC coordinates. (`#4236 <https://github.com/sunpy/sunpy/pull/4236>`__)
- Added `sunpy.coordinates.SunPyBaseCoordinateFrame` and `sunpy.coordinates.BaseHeliographic` to the documentation. (`#4274 <https://github.com/sunpy/sunpy/pull/4274>`__)
- `sunpy.time.TimeRange` had a ``.__contains__`` method and this is now documented. (`#4372 <https://github.com/sunpy/sunpy/pull/4372>`__)
- Revamped sunpy pull request review developer documentation. (`#4378 <https://github.com/sunpy/sunpy/pull/4378>`__)
- Revamped sunpy installation documentation. (`#4378 <https://github.com/sunpy/sunpy/pull/4378>`__)
- Fixed broken documentation links in the guide. (`#4414 <https://github.com/sunpy/sunpy/pull/4414>`__)
- Fixed miscellaneous links in the API documentation. (`#4415 <https://github.com/sunpy/sunpy/pull/4415>`__)
- Added `sunpy.data.data_manager.downloader`, `sunpy.data.data_manager.storage`,
  and `sunpy.net.hek.HEKTable` to the docs. (`#4418 <https://github.com/sunpy/sunpy/pull/4418>`__)
- Added documentation for copying Map objects using the copy module's deepcopy method. (`#4470 <https://github.com/sunpy/sunpy/pull/4470>`__)
- Added information about the :meth:`~sunpy.map.MapSequence.plot` return type. (`#4472 <https://github.com/sunpy/sunpy/pull/4472>`__)
- Added a gallery example for saving and loading sunpy Maps using asdf. (`#4494 <https://github.com/sunpy/sunpy/pull/4494>`__)
- Added description for a counter-intuitive section in the :ref:`sphx_glr_generated_gallery_differential_rotation_reprojected_map.py` example. (`#4548 <https://github.com/sunpy/sunpy/pull/4548>`__)
- Added :ref:`sunpy-topic-guide-coordinates-velocities` to explain how to use velocity information in the coordinates framework. (`#4610 <https://github.com/sunpy/sunpy/pull/4610>`__)
- New gallery example of searching and downloading GOES XRS data (with GOES 15, 16 and 17). (`#4686 <https://github.com/sunpy/sunpy/pull/4686>`__)
- Created the new gallery example :ref:`sphx_glr_generated_gallery_units_and_coordinates_north_offset_frame.py` for `~sunpy.coordinates.NorthOffsetFrame`. (`#4709 <https://github.com/sunpy/sunpy/pull/4709>`__)
- Added more information on which FITS keywords are used for various `sunpy.map.GenericMap`
  properties. (`#4717 <https://github.com/sunpy/sunpy/pull/4717>`__)
- Improved documentation for :func:`sunpy.physics.differential_rotation.diff_rot`. (`#4876 <https://github.com/sunpy/sunpy/pull/4876>`__)


Documentation Fixes
-------------------

- The keyword ``clip_interval`` is now used more extensively in gallery examples when plotting the sample AIA image (e.g., :ref:`sphx_glr_generated_gallery_plotting_aia_example.py`). (`#4573 <https://github.com/sunpy/sunpy/pull/4573>`__)
- Modified :ref:`sphx_glr_generated_gallery_plotting_magnetogram_active_regions.py` to use HMI file from sample data instead of downloading it with Fido. (`#4598 <https://github.com/sunpy/sunpy/pull/4598>`__)
- Removed unnecessary transformations of coordinates prior to plotting them using `~astropy.visualization.wcsaxes.WCSAxes.plot_coord`. (`#4609 <https://github.com/sunpy/sunpy/pull/4609>`__)
- Ensure that all attrs are documented and clean the `sunpy.net.hek.attrs`
  namespace of non-attr objects. (`#4834 <https://github.com/sunpy/sunpy/pull/4834>`__)
- Fixed miscellaneous issues with the gallery example :ref:`sphx_glr_generated_gallery_map_transformations_reprojection_align_aia_hmi.py`. (`#4843 <https://github.com/sunpy/sunpy/pull/4843>`__)
- Fixed the display of arguments in the documentation for `~sunpy.net.Fido` attributes (`sunpy.net.attrs`). (`#4916 <https://github.com/sunpy/sunpy/pull/4916>`__)


Trivial/Internal Changes
------------------------

- ``Fido.fetch`` now always specifies a ``path=`` argument of type `pathlib.Path`
  to the ``fetch`` method of the client. This path will default to the configured
  sunpy download dir, will have the user directory expanded, will have the
  ``{file}`` placeholder and will be tested to ensure that it is writeable. (`#4949 <https://github.com/sunpy/sunpy/pull/4949>`__)
- Added information on what went wrong when `sunpy.map.GenericMap.wcs` fails to parse
  a FITS header into a WCS. (`#4335 <https://github.com/sunpy/sunpy/pull/4335>`__)
- Fixed the `~sunpy.coordinates.frames.Helioprojective` docstring to be clear about the names of the coordinate components. (`#4351 <https://github.com/sunpy/sunpy/pull/4351>`__)
- Raise a better error message if trying to load a FITS file that contains only
  one dimensional data. (`#4426 <https://github.com/sunpy/sunpy/pull/4426>`__)
- The following functions in `sunpy.map` have had their performance greatly increased,
  with runtimes typically improving by a factor of 20x. This has been achieved by
  improving many of the checks so that they only require checking the edge pixels of a
  map as opposed to all of the pixels.

  - :func:`~sunpy.map.contains_full_disk`
  - :func:`~sunpy.map.is_all_off_disk`
  - :func:`~sunpy.map.is_all_on_disk`
  - :func:`~sunpy.map.contains_limb` (`#4463 <https://github.com/sunpy/sunpy/pull/4463>`__)
- Improved the output when you print a sunpy Map. (`#4464 <https://github.com/sunpy/sunpy/pull/4464>`__)
- Creating a `~sunpy.util.MetaDict` with dictionary keys that are not strings now
  raises as user-friendly `ValueError` which prints all the non-compliant keys. (`#4476 <https://github.com/sunpy/sunpy/pull/4476>`__)
- Maps created directly via. `sunpy.map.GenericMap` now have their metadata
  automatically converted to a `~sunpy.util.MetaDict`, which is the same current
  behaviour of the `sunpy.map.Map` factory. (`#4476 <https://github.com/sunpy/sunpy/pull/4476>`__)
- If the ``top_right`` corner given to :meth:`sunpy.map.GenericMap.submap` is
  below or to the right of the ``bottom_left`` corner, a warning is no longer
  raised (as the rectangle is still well defined), but a message is still logged
  at the debug level to the sunpy logger. (`#4491 <https://github.com/sunpy/sunpy/pull/4491>`__)
- Added test support for Python 3.9 (no wheels yet). (`#4569 <https://github.com/sunpy/sunpy/pull/4569>`__)
- ``sunpy.sun`` functions now make use of the `~astropy.coordinates.GeocentricTrueEcliptic` frame to simplify internal calculations, but the returned values are unchanged. (`#4584 <https://github.com/sunpy/sunpy/pull/4584>`__)
- Change the format of the time returned from :func:`~sunpy.coordinates.sun.carrington_rotation_number`
  from ``'jd'`` to ``'iso'``, so printing the `~astropy.time.Time` returned will now print an ISO
  timestamp instead of the Julian days. (`#4819 <https://github.com/sunpy/sunpy/pull/4819>`__)
- The listings for the sample data (`sunpy.data.sample`) are now sorted. (`#4838 <https://github.com/sunpy/sunpy/pull/4838>`__)
- Changed the implementation of a ``hypothesis``-based test so that it does not raise an error with ``hypothesis`` 6.0.0. (`#4852 <https://github.com/sunpy/sunpy/pull/4852>`__)


2.0.0 (2020-06-12)
==================

Backwards Incompatible Changes
------------------------------

- The frames `~sunpy.coordinates.frames.HeliographicStonyhurst` and `~sunpy.coordinates.frames.HeliographicCarrington` now inherit from the new base class `~sunpy.coordinates.frames.BaseHeliographic`.
  This changes means that ``isinstance(frame, HeliographicStonyhurst)`` is no longer ``True`` when ``frame`` is `~sunpy.coordinates.frames.HeliographicCarrington`. (`#3595 <https://github.com/sunpy/sunpy/pull/3595>`__)
- `~sunpy.visualization.colormaps.color_tables.aia_color_table`, `~sunpy.visualization.colormaps.color_tables.eit_color_table` and `~sunpy.visualization.colormaps.color_tables.suvi_color_table` now only take `astropy.units` quantities instead of strings. (`#3640 <https://github.com/sunpy/sunpy/pull/3640>`__)
- `sunpy.map.Map` is now more strict when the metadata of a map cannot be validated, and
  an error is now thrown instead of a warning if metadata cannot be validated. In order to
  load maps that previously loaded without error you may need to pass ``silence_errors=True``
  to `sunpy.map.Map`. (`#3646 <https://github.com/sunpy/sunpy/pull/3646>`__)
- ``Fido.search`` will now return results from all clients which match a query, you no longer have to make the query specific to a single client. This means that searches involving the 'eve' and 'rhessi' instruments will now potentially also return results from the VSO. For `~sunpy.net.dataretriever.RHESSIClient` you can now specify ``a.Physobs("summary_lightcurve")`` to only include the summary lightcurve data products not provided by the VSO. (`#3770 <https://github.com/sunpy/sunpy/pull/3770>`__)
- The objects returned by the ``search`` methods on ``VSOClient``, ``JSOCClient`` and ``GenericClient`` have been changed to be based on `sunpy.net.base_client.BaseQueryResponse`. This introduces a few subtle breaking changes for people using the client search methods directly (not ``Fido.search``), or people using ``sunpy.net.fido_factory.UnifiedResponse.get_response``. When slicing an instance of ``QueryResponse`` it will now return an instance of itself, ``QueryResponse.blocks`` can be used to access the underlying records. Also, the ``.client`` attribute of the response no longer has to be the instance of the class the search was made with, however, it often is. (`#3770 <https://github.com/sunpy/sunpy/pull/3770>`__)
- `~sunpy.coordinates.frames.HeliographicCarrington` is now an observer-based frame, where the observer location (specifically, the distance from the Sun) is used to account for light travel time when determining apparent Carrington longitudes.  Coordinate transformations using this frame now require an observer to be specified. (`#3782 <https://github.com/sunpy/sunpy/pull/3782>`__)
- To enable the precise co-alignment of solar images from different observatories, the calculation of Carrington coordinates now ignores the stellar-aberration correction due to observer motion.
  For an Earth observer, this change results in an increase in Carrington longitude of ~20 arcseconds.
  See :ref:`sunpy-topic-guide-coordinates-carrington` for more information. (`#3782 <https://github.com/sunpy/sunpy/pull/3782>`__)
- Fixed a bug where some of the coordinate transformations could raise `ValueError` instead of `~astropy.coordinates.ConvertError` when the transformation could not be performed. (`#3894 <https://github.com/sunpy/sunpy/pull/3894>`__)
- Astropy 3.2 is now the minimum required version of that dependency. (`#3936 <https://github.com/sunpy/sunpy/pull/3936>`__)


Deprecations and Removals
-------------------------

- Fido search attrs available as `sunpy.net.attrs` i.e, ``a.Time``, ``a.Instrument`` etc are now deprecated as VSO attrs (`sunpy.net.vso.attrs`). (`#3714 <https://github.com/sunpy/sunpy/pull/3714>`__)
- ``sunpy.util.multimethod.MultiMethod`` is deprecated, `functools.singledispatch` provides equivalent functionality in the standard library. (`#3714 <https://github.com/sunpy/sunpy/pull/3714>`__)
- ``sunpy.net.vso.attrs.Physobs`` has been moved to `sunpy.net.attrs.Physobs` and the original deprecated. (`#3877 <https://github.com/sunpy/sunpy/pull/3877>`__)
- Deprecate ``sunpy.instr.aia.aiaprep`` in favor of the `aiapy.calibrate.register` function in the
  [aiapy](https://gitlab.com/LMSAL_HUB/aia_hub/aiapy) package.
  ``sunpy.instr.aia.aiaprep`` will be removed in version 2.1 (`#3960 <https://github.com/sunpy/sunpy/pull/3960>`__)
- Removed the module ``sunpy.sun.sun``, which was deprecated in version 1.0.
  Use the module `sunpy.coordinates.sun` instead. (`#4014 <https://github.com/sunpy/sunpy/pull/4014>`__)
- Removed Sun-associated functions in `sunpy.coordinates.ephemeris`, which were deprecated in 1.0.
  Use the corresponding functions in `sunpy.coordinates.sun`. (`#4014 <https://github.com/sunpy/sunpy/pull/4014>`__)
- Remove the deprecated ``sunpy.net.vso.vso.VSOClient`` ``.query_legacy`` and ``.latest`` methods. (`#4109 <https://github.com/sunpy/sunpy/pull/4109>`__)
- Removed the sample datasets NOAAINDICES_TIMESERIES and NOAAPREDICT_TIMESERIES because they will invariably be out of date.
  Up-to-date versions of these NOAA indices can be downloaded using `~sunpy.net.Fido` (see :ref:`sphx_glr_generated_gallery_plotting_solar_cycle_example.py`). (`#4169 <https://github.com/sunpy/sunpy/pull/4169>`__)


Features
--------

- Added `~sunpy.coordinates.metaframes.RotatedSunFrame` for defining coordinate frames that account for solar rotation. (`#3537 <https://github.com/sunpy/sunpy/pull/3537>`__)
- Added a context manager (`~sunpy.coordinates.transform_with_sun_center`) to ignore the motion of the center of the Sun for coordinate transformations. (`#3540 <https://github.com/sunpy/sunpy/pull/3540>`__)
- Updated the gallery example titled 'Downloading and plotting an HMI magnetogram' to rotate the HMI magnetogram such that solar North is pointed up. (`#3573 <https://github.com/sunpy/sunpy/pull/3573>`__)
- Creates a function named ``sunpy.map.sample_at_coords`` that samples the data from the map at the given set of coordinates. (`#3592 <https://github.com/sunpy/sunpy/pull/3592>`__)
- Enabled the discovery of search attributes for each of our clients. (`#3637 <https://github.com/sunpy/sunpy/pull/3637>`__)
- Printing `sunpy.net.attrs.Instrument` or other "attrs" will show all attributes that exist under the corresponding "attr". (`#3637 <https://github.com/sunpy/sunpy/pull/3637>`__)
- Printing `sunpy.net.Fido` will print out all the clients that Fido can use. (`#3637 <https://github.com/sunpy/sunpy/pull/3637>`__)
- Updates `~sunpy.map.GenericMap.draw_grid` to allow disabling the axes labels and the ticks on the top and right axes. (`#3673 <https://github.com/sunpy/sunpy/pull/3673>`__)
- Creates a ``tables`` property for `~sunpy.net.fido_factory.UnifiedResponse`, which allows to access the `~sunpy.net.base_client.BaseQueryResponse` as an `~astropy.table.Table`, which then can be used for indexing of results. (`#3675 <https://github.com/sunpy/sunpy/pull/3675>`__)
- Change the APIs for ``sunpy.map.GenericMap.draw_rectangle`` and :meth:`sunpy.map.GenericMap.submap` to be consistent with each other and to use keyword-only arguments for specifying the bounding box. (`#3677 <https://github.com/sunpy/sunpy/pull/3677>`__)
- Updates the `~sunpy.map.GenericMap.observer_coordinate` property to warn the user of specific missing metadata for each frame.
  Omits warning about frames where all metadata is missing or all meta is present. (`#3692 <https://github.com/sunpy/sunpy/pull/3692>`__)
- Added `sunpy.util.config.copy_default_config` that copies the default config file to the user's config directory. (`#3722 <https://github.com/sunpy/sunpy/pull/3722>`__)
- ``sunpy.database`` now supports adding database entries and downloading data from ``HEK`` query (`#3731 <https://github.com/sunpy/sunpy/pull/3731>`__)
- Added a helper function (`~sunpy.coordinates.utils.get_rectangle_coordinates`) for defining a rectangle in longitude and latitude coordinates. (`#3737 <https://github.com/sunpy/sunpy/pull/3737>`__)
- Add a ``.data`` property in `~sunpy.timeseries.GenericTimeSeries`, so that users are encouraged to use :meth:`~sunpy.timeseries.GenericTimeSeries.to_dataframe` to get the data of the timeseries. (`#3746 <https://github.com/sunpy/sunpy/pull/3746>`__)
- It is now possible to turn on or off various corrections in :func:`~sunpy.coordinates.sun.L0` (the apparent Carrington longitude of Sun-disk center as seen from Earth). (`#3782 <https://github.com/sunpy/sunpy/pull/3782>`__)
- Made skimage.transform import lazy to reduce import time of `sunpy.image.transform` by ~50% (`#3818 <https://github.com/sunpy/sunpy/pull/3818>`__)
- Add support for parfive 1.1. This sets a limit on the number of open connections to JSOC when downloading files to 10. (`#3822 <https://github.com/sunpy/sunpy/pull/3822>`__)
- Fido clients (subclasses of `sunpy.net.base_client.BaseClient`) can now register their own attrs modules with `sunpy.net.attrs`.
  This allows clients which require attr classes specific to that client to register modules that can be used by the user i.e. ``a.vso``.
  It also allows clients implemented externally to sunpy to register attrs. (`#3869 <https://github.com/sunpy/sunpy/pull/3869>`__)
- Added the methods :meth:`sunpy.map.GenericMap.quicklook` and :meth:`sunpy.map.MapSequence.quicklook` to display an HTML summary of the instance, including interactive controls.
  When using Jupyter notebooks, this HTML summary is automatically shown instead of a text-only representation. (`#3951 <https://github.com/sunpy/sunpy/pull/3951>`__)
- Added `_localfilelist` method in ``sunpy.util.scraper.Scraper`` to scrap local data archives. (`#3994 <https://github.com/sunpy/sunpy/pull/3994>`__)
- Added extra constants to `sunpy.sun.constants`:

  - Longitude of the prime meridian (epoch J2000.0) : ``sunpy.sun.constants.get('W_0')``
  - Sidereal rotation rate : `sunpy.sun.constants.sidereal_rotation_rate`
  - First Carrington rotation (JD TT) : `sunpy.sun.constants.first_carrington_rotation`
  - Mean synodic period : `sunpy.sun.constants.mean_synodic_period`
  - Right ascension (RA) of the north pole (epoch J2000.0) : ``sunpy.sun.constants.get('alpha_0')``
  - Declination of the north pole (epoch J2000.0) : ``sunpy.sun.constants.get('delta_0')`` (`#4013 <https://github.com/sunpy/sunpy/pull/4013>`__)
- Adds to ``sunpy.util.scraper.Scraper`` the ability to include regular expressions in the URL passed. (`#4107 <https://github.com/sunpy/sunpy/pull/4107>`__)


Bug Fixes
---------

- Added support for passing ``TimeSeriesMetaData`` object to ``timeseries_factory`` and associated validation tests. (`#3639 <https://github.com/sunpy/sunpy/pull/3639>`__)
- Now when `~sunpy.map.GenericMap` fails to load a file, the filename that failed to load will now be part of the error message. (`#3727 <https://github.com/sunpy/sunpy/pull/3727>`__)
- Work around incorrect Content-Disposition headers in some VSO downloads, which were leading to mangled filenames. (`#3740 <https://github.com/sunpy/sunpy/pull/3740>`__)
- ``Fido.search`` can now service queries without ``a.Time`` being specified. This is currently only used by the `sunpy.net.jsoc.JSOCClient`. (`#3770 <https://github.com/sunpy/sunpy/pull/3770>`__)
- Fixed a bug with the calculation of Carrington longitude as seen from Earth where it was using an old approach instead of the current approach (for example, the varying Sun-Earth distance is now taken into account).
  The old approach resulted in errors no greater than 7 arcseconds in Carrington longitude when using `~sunpy.coordinates.sun.L0` and `~sunpy.coordinates.frames.HeliographicCarrington`. (`#3772 <https://github.com/sunpy/sunpy/pull/3772>`__)
- Updated `sunpy.map.CompositeMap.plot` to support a linewidths argument. (`#3792 <https://github.com/sunpy/sunpy/pull/3792>`__)
- Fix a bug in `sunpy.net.jsoc.JSOCClient` where requesting data for export would not work if a non-time primekey was used. (`#3825 <https://github.com/sunpy/sunpy/pull/3825>`__)
- Add support for passing paths of type `pathlib.Path` in `sunpy.net.jsoc.JSOCClient.fetch`. (`#3838 <https://github.com/sunpy/sunpy/pull/3838>`__)
- Add explicit support for dealing with download urls for files, under 'as-is' protocol in `sunpy.net.jsoc.JSOCClient.get_request`. (`#3838 <https://github.com/sunpy/sunpy/pull/3838>`__)
- Updated the method used to filter time in the VSO post-search filtering function. (`#3840 <https://github.com/sunpy/sunpy/pull/3840>`__)
- Fix failing of fetching of the indexed JSOCResponses using `sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch`. (`#3852 <https://github.com/sunpy/sunpy/pull/3852>`__)
- Prevented `sunpy.map.GenericMap.plot` modifying in-place any items passed as ``imshow_kwargs``. (`#3867 <https://github.com/sunpy/sunpy/pull/3867>`__)
- Changed the format of DATE-OBS in `sunpy.map.GenericMap.wcs` from iso to isot (ie. with a "T" between the date and time) to conform with the FITS standard. (`#3872 <https://github.com/sunpy/sunpy/pull/3872>`__)
- Fixed a minor error (up to ~10 arcseconds) in the calculation of the Sun's position angle (:func:`sunpy.coordinates.sun.P`). (`#3886 <https://github.com/sunpy/sunpy/pull/3886>`__)
- `~sunpy.net.hek.HEKClient` was returning HTML and not JSON. (`#3899 <https://github.com/sunpy/sunpy/pull/3899>`__)
- Updated to HTTPS for HEK. (`#3917 <https://github.com/sunpy/sunpy/pull/3917>`__)
- The accuracy of the output of :func:`sunpy.coordinates.ephemeris.get_horizons_coord` is significantly improved. (`#3919 <https://github.com/sunpy/sunpy/pull/3919>`__)
- Fixed a bug where the longitude value for the reference coordinate in the Map repr would be displayed with the unintended longitude wrapping. (`#3959 <https://github.com/sunpy/sunpy/pull/3959>`__)
- It is now possible to specify a local file path to
  `sunpy.data.data_manager.DataManager.override_file` without having to prefix it
  with ``file://``. (`#3970 <https://github.com/sunpy/sunpy/pull/3970>`__)
- Closed the session in the destructor of VSOClient thus solving the problem of socket being left open (`#3973 <https://github.com/sunpy/sunpy/pull/3973>`__)
- Fixed a bug of where results of VSO searches would have inconsistent ordering in ``sunpy.net.vso.vso.QueryResponse`` by always sorting the results by start time. (`#3974 <https://github.com/sunpy/sunpy/pull/3974>`__)
- Fixes two bugs in `sunpy.util.deprecated`: correctly calculates the
  removal version and does not override the default and/or alternative functionality
  message. Providing a custom deprecation message now suppresses any
  mention of the removal version. Additionally, a ``pending`` keyword argument is
  provided to denote functions/classes that are pending deprecation. (`#3982 <https://github.com/sunpy/sunpy/pull/3982>`__)
- Correctly generate labels for sliders in
  ``~sunpy.visualization.animator.ArrayAnimatorWCS`` when the number of pixel
  dimensions and the number of world dimensions are not the same in the WCS. (`#3990 <https://github.com/sunpy/sunpy/pull/3990>`__)
- Updated VSOClient.response_block_properties to check if "None" is in the return. (`#3993 <https://github.com/sunpy/sunpy/pull/3993>`__)
- Fix a bug with ``sunpy.visualization.animator.ArrayAnimatorWCS`` where animating
  a line with a masked array with the whole of the initial line masked out the
  axes limits for the x axis were not correctly set. (`#4001 <https://github.com/sunpy/sunpy/pull/4001>`__)
- Fixed passing in a list of URLs into `sunpy.map.GenericMap`, before it caused an error due to the wrong type being returned. (`#4007 <https://github.com/sunpy/sunpy/pull/4007>`__)
- Fixed a bug with :func:`~sunpy.coordinates.transform_with_sun_center` where the global variable was sometimes restored incorrectly.
  This bug was most likely encountered if there was a nested use of this context manager. (`#4015 <https://github.com/sunpy/sunpy/pull/4015>`__)
- Fixes a bug in fido_factory to allow  path="./" in fido.fetch(). (`#4058 <https://github.com/sunpy/sunpy/pull/4058>`__)
- Prevented ``sunpy.io.fits.header_to_fits`` modifying the passed header in-place. (`#4067 <https://github.com/sunpy/sunpy/pull/4067>`__)
- Strip out any unknown unicode from the HEK response to prevent it failing to load some results. (`#4088 <https://github.com/sunpy/sunpy/pull/4088>`__)
- Fixed a bug in :func:`~sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst` that resulted in a error when requesting an array of locations in conjunction with enabling the light-travel-time correction. (`#4112 <https://github.com/sunpy/sunpy/pull/4112>`__)
- `sunpy.map.GenericMap.top_right_coord` and `~sunpy.map.GenericMap.center`
  have had their definitions clarified, and both have had off-by-one indexing
  errors fixed. (`#4121 <https://github.com/sunpy/sunpy/pull/4121>`__)
- Fixed `sunpy.map.GenericMap.submap()` when scaled pixel units (e.g. ``u.mpix``)
  are used. (`#4127 <https://github.com/sunpy/sunpy/pull/4127>`__)
- Fixed bugs in ``sunpy.util.scraper.Scraper.filelist``
  that resulted in error when the HTML page of URL opened by the scraper contains some "a" tags without "href" attribute
  and resulted in incorrect file urls when any href stores filepath relative to the URL's domain instead of just a filename. (`#4132 <https://github.com/sunpy/sunpy/pull/4132>`__)
- Fixed inconsistencies in how `~sunpy.map.GenericMap.submap` behaves when passed corners in pixel and world coordinates.
  The behavior for submaps specified in pixel coordinates is now well-defined for pixels on the boundary of the rectangle
  and is consistent for all boundaries. Previously pixels on the lower left boundary were included, but excluded on the
  upper and right boundary. This means the shape of a submap may now be 1 pixel larger in each dimension.
  Added several more tests for `~sunpy.map.GenericMap.submap` for a range of cutout sizes in both pixel and world
  coordinates. (`#4134 <https://github.com/sunpy/sunpy/pull/4134>`__)
- `sunpy.map.on_disk_bounding_coordinates` now fully propagates the coordinate
  frame of the input map to the output coordinates. Previously only the observer
  coordinate, and no other frame attributes, were propagated. (`#4141 <https://github.com/sunpy/sunpy/pull/4141>`__)
- Fix an off-by-one error in the reference pixel returned by
  :func:`sunpy.map.header_helper.make_fitswcs_header`. (`#4152 <https://github.com/sunpy/sunpy/pull/4152>`__)
- `sunpy.map.GenericMap.reference_pixel` now uses zero-based indexing, in order
  to be consistent with the rest of the `sunpy.map` API. (`#4154 <https://github.com/sunpy/sunpy/pull/4154>`__)
- Previously `sunpy.map.GenericMap.resample` with ``method='linear'`` was
  using an incorrect and constant value to fill edges when upsampling a map. Values
  near the edges are now correctly extrapolated using the ``fill_value=extrapolate``
  option to `scipy.interpolate.interp1d`. (`#4164 <https://github.com/sunpy/sunpy/pull/4164>`__)
- Fixed a bug where passing an `int` or `list` via the ``hdus`` keyword argument to
  ``sunpy.io.fits.read`` threw an exception because the list of HDU objects was no longer
  of type `~astropy.io.fits.HDUList`. (`#4183 <https://github.com/sunpy/sunpy/pull/4183>`__)
- Fix attr printing when the attr registry is empty for that attr (`#4199 <https://github.com/sunpy/sunpy/pull/4199>`__)
- Improved the accuracy of :func:`~sunpy.coordinates.sun.angular_radius` by removing the use of the small-angle approximation.
  The inaccuracy had been less than 5 milliarcseconds. (`#4239 <https://github.com/sunpy/sunpy/pull/4239>`__)
- Fixed a bug with the ``observer`` frame attribute for coordinate frames where an input that was not supplied as a `~astropy.coordinates.SkyCoord` would sometimes result in a transformation error. (`#4266 <https://github.com/sunpy/sunpy/pull/4266>`__)


Improved Documentation
----------------------

- Fixed an issue with the scaling of class-inheritance diagrams in the online documentation by blocking the versions of graphviz containing a bug. (`#3548 <https://github.com/sunpy/sunpy/pull/3548>`__)
- A new example gallery example "Plotting a difference image" has been added,
  which can be used for base difference or running difference images. (`#3627 <https://github.com/sunpy/sunpy/pull/3627>`__)
- Removed obsolete Astropy Helpers submodule section in :file:`CONTRIBUTING.rst`;
  Also removed mentions of astropy_helpers in all files of the project. (`#3676 <https://github.com/sunpy/sunpy/pull/3676>`__)
- Corrected misleading `~sunpy.timeseries.metadata.TimeSeriesMetaData` documentation about optional parameters. (`#3680 <https://github.com/sunpy/sunpy/pull/3680>`__)
- Added an example for `~sunpy.map.GenericMap.world_to_pixel` function in the Units & Coordinates guide. (`#3776 <https://github.com/sunpy/sunpy/pull/3776>`__)
- Added a :ref:`page <sunpy-topic-guide-coordinates-carrington>` describing how SunPy calculates Carrington longitudes. (`#3782 <https://github.com/sunpy/sunpy/pull/3782>`__)
- Changed padding value of an example in the example gallery to fix the overlap of titles and x-label axes. (`#3835 <https://github.com/sunpy/sunpy/pull/3835>`__)
- More information and links about how to create changelogs. (`#3856 <https://github.com/sunpy/sunpy/pull/3856>`__)
- Clarified some inputs to `sunpy.map.GenericMap.plot`. (`#3866 <https://github.com/sunpy/sunpy/pull/3866>`__)
- Changed quoted sentence (that we suggest authors add to their research papers) in CITATION.rst (`#3896 <https://github.com/sunpy/sunpy/pull/3896>`__)
- Add example of how to use SunPy's HEK client to search for the GOES flare event list. (`#3953 <https://github.com/sunpy/sunpy/pull/3953>`__)
- Improved the doc layout of `sunpy.data.sample`. (`#4034 <https://github.com/sunpy/sunpy/pull/4034>`__)
- Made improvements to STEREO starfield gallery example. (`#4039 <https://github.com/sunpy/sunpy/pull/4039>`__)
- Improved the documentation of `sunpy.map.GenericMap.resample`. (`#4043 <https://github.com/sunpy/sunpy/pull/4043>`__)
- Updated the STEREO starfield example to use all of the information in the star catalog. (`#4116 <https://github.com/sunpy/sunpy/pull/4116>`__)
- Mini-galleries are now easier to create in the documentation thanks to a custom Sphinx directive (``minigallery``).
  The page :ref:`sunpy-topic-guide-coordinates-rotatedsunframe` has an example of a mini-gallery at the bottom. (`#4124 <https://github.com/sunpy/sunpy/pull/4124>`__)
- Added `sunpy.visualization.colormaps.color_tables` to the docs. (`#4182 <https://github.com/sunpy/sunpy/pull/4182>`__)
- Made minor improvements to the map histogramming example. (`#4205 <https://github.com/sunpy/sunpy/pull/4205>`__)
- Add a warning to ``sunpy.io`` docs to recommend not using it for FITS (`#4208 <https://github.com/sunpy/sunpy/pull/4208>`__)


Trivial/Internal Changes
------------------------

- Removed un-used and un-tested code paths in the private ``_remove_lytaf_events`` function
  in ``sunpy.instr.lyra``. (`#3570 <https://github.com/sunpy/sunpy/pull/3570>`__)
- Removed ``astropy_helpers`` and this means that ``python setup.py <test,build_docs>`` no longer works.
  So if you want to:

  * Run the tests: Use ``tox -e <env name>`` or call ``pytest`` directly
  * Build the docs: Use ``tox -e docs`` or cd into the docs folder and run ``make html`` or ``sphinx-build docs docs/_build/html -W -b html -d docs/_build/.doctrees`` (`#3598 <https://github.com/sunpy/sunpy/pull/3598>`__)
- Cleaned up test warnings in sunpy.coordinates. (`#3652 <https://github.com/sunpy/sunpy/pull/3652>`__)
- Fix Python version for requiring importlib_resources (`#3683 <https://github.com/sunpy/sunpy/pull/3683>`__)
- `sunpy.net.attr.AttrWalker` no longer uses ``sunpy.util.multimethod.MultiMethod`` it uses a derivative of `functools.singledispatch` `sunpy.util.functools.seconddispatch` which dispatches on the second argument. (`#3714 <https://github.com/sunpy/sunpy/pull/3714>`__)
- Errors from a VSO search will now be raised to the user. (`#3719 <https://github.com/sunpy/sunpy/pull/3719>`__)
- Fixed the transformation test for `~sunpy.coordinates.metaframes.NorthOffsetFrame`, which would intermittently fail. (`#3775 <https://github.com/sunpy/sunpy/pull/3775>`__)
- :func:`~sunpy.coordinates.sun.earth_distance` is now computed without using coordinate transformations for better performance. (`#3782 <https://github.com/sunpy/sunpy/pull/3782>`__)
- Created a helper function for testing the equality/closeness of longitude angles (i.e., angles with wrapping). (`#3804 <https://github.com/sunpy/sunpy/pull/3804>`__)
- Bump the astropy version figure tests are run with from 3.1.2 to 3.2.3 (`#3925 <https://github.com/sunpy/sunpy/pull/3925>`__)
- Used `urllib.parse.urlsplit` in ``sunpy.util.scraper`` for file scraping functionalities. (`#3956 <https://github.com/sunpy/sunpy/pull/3956>`__)
- Added `sunpy.net.base_client.BaseClient.check_attr_types_in_query` as a helper method
  to check if a query contains a set of required attributes, and is a subset of optional
  attributes. (`#3979 <https://github.com/sunpy/sunpy/pull/3979>`__)
- Removes appending login details for ftp urls from scraper. (`#4020 <https://github.com/sunpy/sunpy/pull/4020>`__)
- Re-factored the `sunpy.map.Map` factory to dispatch argument parsing based on type. (`#4037 <https://github.com/sunpy/sunpy/pull/4037>`__)
- Improved the error message raised by the Map factory when a map matches multiple source map types. (`#4052 <https://github.com/sunpy/sunpy/pull/4052>`__)
- Added log messages when the sample data fails to download. (`#4137 <https://github.com/sunpy/sunpy/pull/4137>`__)
- Remove an Astropy 3.1 compatibility wrapper for ``Quantity.to_string``. (`#4172 <https://github.com/sunpy/sunpy/pull/4172>`__)
- Refactor the sphinx config to no longer depend on astropy-sphinx and more
  closely match the new sunpy package template (`#4188 <https://github.com/sunpy/sunpy/pull/4188>`__)


1.1.0 (2020-01-10)
==================

Backwards Incompatible Changes
------------------------------

- The ``sunpy.net.vso.vso.get_online_vso_url`` function has been broken into two components, the new ``sunpy.net.vso.vso.get_online_vso_url`` function takes no arguments (it used to take three) and now only returns an online VSO mirror or None.
  The construction of a ``zeep.Client`` object is now handled by ``sunpy.net.vso.vso.build_client`` which has a more flexible API for customising the ``zeep.Client`` interface. (`#3330 <https://github.com/sunpy/sunpy/pull/3330>`__)
- Importing ``sunpy.timeseries.timeseriesbase`` no longer automatically imports
  Matplotlib. (`#3376 <https://github.com/sunpy/sunpy/pull/3376>`__)
- :meth:`sunpy.timeseries.sources.NOAAIndicesTimeSeries.peek()` now checks that the `type` argument is a
  valid string, and raises a `ValueError` if it isn't. (`#3378 <https://github.com/sunpy/sunpy/pull/3378>`__)
- Observer-based coordinate frames (`~sunpy.coordinates.frames.Heliocentric` and `~sunpy.coordinates.frames.Helioprojective`) no longer assume a default observer (Earth) if no observer is specified.  These frames can now be used with no observer specified, but most transformations cannot be performed for such frames.  This removal of a default observer only affects `sunpy.coordinates`, and has no impact on the default observer in `sunpy.map`. (`#3388 <https://github.com/sunpy/sunpy/pull/3388>`__)
- The callback functions provided to
  ``BaseFuncAnimator`` ``button_func`` keyword
  argument now take two positional arguments rather than one. The function
  signature is now ``(animator, event)`` where the first arg is the animator
  object, and the second is the matplotlib mouse event. (`#3407 <https://github.com/sunpy/sunpy/pull/3407>`__)
- The colormap stored in SunPy's Map subclasses (ie. ``map.plot_settings['cmap']``)
  can now be colormap string instead of the full `matplotlib.colors.Colormap`
  object. To get the full `~matplotlib.colors.Colormap` object use the new attribute
  ``map.cmap``. (`#3412 <https://github.com/sunpy/sunpy/pull/3412>`__)
- Fix a warning in `sunpy.map.GenericMap.rotate` where the truth value of an array
  was being calculated. This changes the behaviour of
  `~sunpy.map.GenericMap.rotate` when the ``angle=`` parameter is not an
  `~astropy.units.Quantity` object to raise `TypeError` rather than `ValueError`. (`#3456 <https://github.com/sunpy/sunpy/pull/3456>`__)


Deprecations and Removals
-------------------------

- Removed the step of reparing images (replacing non-finite entries with local mean) before coaligning them. The user is expected to do this themselves before coaligning images. If NaNs/non-finite entries are present, a warning is thrown.
  The function ``sunpy.image.coalignment.repair_image_nonfinite`` is deprecated. (`#3287 <https://github.com/sunpy/sunpy/pull/3287>`__)
- The method to convert a `~sunpy.coordinates.frames.Helioprojective` frame from 2D to 3D has been renamed from ``calculate_distance`` to `~sunpy.coordinates.frames.Helioprojective.make_3d`.  This method is not typically directly called by users. (`#3389 <https://github.com/sunpy/sunpy/pull/3389>`__)
- ``sunpy.visualization.animator.ImageAnimatorWCS`` is now deprecated in favour of
  ``ArrayAnimatorWCS``. (`#3407 <https://github.com/sunpy/sunpy/pull/3407>`__)
- ``sunpy.cm`` has been moved to `sunpy.visualization.colormaps` and will be
  removed in a future version. (`#3410 <https://github.com/sunpy/sunpy/pull/3410>`__)


Features
--------

- Add a new `sunpy.data.manager` and `sunpy.data.cache` for dealing with versioned remote data within functions.
  Please see the ``Remote Data Manager`` guide. (`#3124 <https://github.com/sunpy/sunpy/pull/3124>`__)
- Added the coordinate frames `~sunpy.coordinates.frames.HeliocentricEarthEcliptic` (HEE), `~sunpy.coordinates.frames.GeocentricSolarEcliptic` (GSE), `~sunpy.coordinates.frames.HeliocentricInertial` (HCI), and `~sunpy.coordinates.frames.GeocentricEarthEquatorial` (GEI). (`#3212 <https://github.com/sunpy/sunpy/pull/3212>`__)
- Added SunPy Map support for GOES SUVI images. (`#3269 <https://github.com/sunpy/sunpy/pull/3269>`__)
- - Support APE14 for ``ImageAnimatorWCS`` in SunPy's visualization module (`#3275 <https://github.com/sunpy/sunpy/pull/3275>`__)
- Add ability to disable progressbars when downloading files using ``sunpy.net.helioviewer`` and edited docstrings to mention this feature. (`#3280 <https://github.com/sunpy/sunpy/pull/3280>`__)
- Adds support for searching and downloading SUVI data. (`#3301 <https://github.com/sunpy/sunpy/pull/3301>`__)
- Log all VSO XML requests and responses to the SunPy logger at the ``DEBUG``
  level. (`#3330 <https://github.com/sunpy/sunpy/pull/3330>`__)
- Transformations between frames in `sunpy.coordinates` can now provide detailed debugging output.  Set the `logging` level to ``DEBUG`` to enable this output. (`#3339 <https://github.com/sunpy/sunpy/pull/3339>`__)
- Added the `sunpy.coordinates.sun.carrington_rotation_time` function to
  compute the time of a given Carrington rotation number. (`#3360 <https://github.com/sunpy/sunpy/pull/3360>`__)
- A new method has been added to remove columns from a
  `sunpy.timeseries.GenericTimeSeries`. (`#3361 <https://github.com/sunpy/sunpy/pull/3361>`__)
- Add ``shape`` property to TimeSeries. (`#3380 <https://github.com/sunpy/sunpy/pull/3380>`__)
- Added ASDF schemas for the new coordinate frames (`~sunpy.coordinates.frames.GeocentricEarthEquatorial`, `~sunpy.coordinates.frames.GeocentricSolarEcliptic`, `~sunpy.coordinates.frames.HeliocentricEarthEcliptic`, `~sunpy.coordinates.frames.HeliocentricInertial`).  See the gallery for an example of using ``asdf`` to save and load a coordinate frame. (`#3398 <https://github.com/sunpy/sunpy/pull/3398>`__)
- ``sunpy.visualization.animator.ArrayAnimatorWCS`` was added which uses the WCS
  object to get the coordinates of all axes, including the slider labels. It also provides the
  ability to customise the plot by specifying arguments to
  `~astropy.visualization.wcsaxes.WCSAxes` methods and supports animation of
  WCS aware line plots with Astroy 4.0. (`#3407 <https://github.com/sunpy/sunpy/pull/3407>`__)
- The returned list of `~sunpy.map.Map` objects is now sorted by filename when
  passing a directory or glob pattern to `~sunpy.map.map_factory.MapFactory`. (`#3408 <https://github.com/sunpy/sunpy/pull/3408>`__)
- Single character wildcards and character ranges can now be passed as
  glob patterns to `~sunpy.map.Map`. (`#3408 <https://github.com/sunpy/sunpy/pull/3408>`__)
- `~sunpy.map.Map` now accepts filenames and directories as `pathlib.Path`
  objects. (`#3408 <https://github.com/sunpy/sunpy/pull/3408>`__)
- `~sunpy.map.GenericMap` objects now have a ``.cmap`` attribute, which returns the full `~matplotlib.colors.Colormap`.
  object. (`#3412 <https://github.com/sunpy/sunpy/pull/3412>`__)
- ``sunpy.io.write_file`` now accepts `~pathlib.Path` objects as filename inputs. (`#3469 <https://github.com/sunpy/sunpy/pull/3469>`__)
- :func:`sunpy.map.header_helper.make_fitswcs_header` now accepts a `tuple` representing the shape of an array as well as the actual array as the ``data`` argument. (`#3483 <https://github.com/sunpy/sunpy/pull/3483>`__)
- Made a couple of module imports lazy to reduce the import time of sunpy.map by
  ~40%. (`#3495 <https://github.com/sunpy/sunpy/pull/3495>`__)
- `sunpy.map.GenericMap.wcs` now uses the full FITS header to construct the WCS.
  This adds support for instruments with more complex projections, such as WISPR,
  however does mean that Map will be more sensitive to incorrect or invalid FITS
  headers. If you are using custom headers with SunPy Map you might encounter
  issues relating to this change. (`#3501 <https://github.com/sunpy/sunpy/pull/3501>`__)
- ``sunpy.visualization.animator.BaseFuncAnimator`` now takes an optional
  ``slider_labels`` keyword argument which draws text labels in the center of the
  sliders. (`#3504 <https://github.com/sunpy/sunpy/pull/3504>`__)
- Added a more helpful error message when trying to load a file or directory
  that doesn't exist with `sunpy.map.Map`. (`#3568 <https://github.com/sunpy/sunpy/pull/3568>`__)
- Add ``__repr__`` for `~sunpy.map.MapSequence` objects  so that users can view the
  critical information of all the ``Map`` objects, in a concise manner. (`#3636 <https://github.com/sunpy/sunpy/pull/3636>`__)


Bug Fixes
---------

- Fixed accuracy issues with the calculations of Carrington longitude (`~sunpy.coordinates.sun.L0`) and Carrington rotation number (`~sunpy.coordinates.sun.carrington_rotation_number`). (`#3178 <https://github.com/sunpy/sunpy/pull/3178>`__)
- Updated :func:`sunpy.map.header_helper.make_fitswcs_header` to be more strict on the inputs it accepts. (`#3183 <https://github.com/sunpy/sunpy/pull/3183>`__)
- Fix the calculation of ``rsun_ref`` in :func:`~sunpy.map.header_helper.make_fitswcs_header` and and
  ensure that the default reference pixel is indexed from 1. (`#3184 <https://github.com/sunpy/sunpy/pull/3184>`__)
- Fixed the missing transformation between two `~sunpy.coordinates.HeliographicCarrington` frames with different observation times. (`#3186 <https://github.com/sunpy/sunpy/pull/3186>`__)
- `sunpy.map.sources.AIAMap` and `sunpy.map.sources.HMIMap` will no longer assume
  the existence of certain header keys. (`#3217 <https://github.com/sunpy/sunpy/pull/3217>`__)
- :func:`sunpy.map.header_helper.make_fitswcs_header` now supports specifying the map projection
  rather than defaulting to ``TAN``. (`#3218 <https://github.com/sunpy/sunpy/pull/3218>`__)
- Fix the behaviour of
  ``sunpy.coordinates.frames.Helioprojective.calculate_distance`` if the
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
- Fixed all coordinate transformations to properly handle a change in observation time. (`#3247 <https://github.com/sunpy/sunpy/pull/3247>`__)
- Fixed the handling of coordinates with velocity information when transforming between Astropy frames and SunPy frames. (`#3247 <https://github.com/sunpy/sunpy/pull/3247>`__)
- Fixed ``sunpy.physics.solar_rotation.calculate_solar_rotate_shift`` so that it does not calculate a shift between the reference layer and itself, which would sometimes incorrectly result in a shift of a pixel due to numerical precision. (`#3255 <https://github.com/sunpy/sunpy/pull/3255>`__)
- Stop crash when ``LineAnimator`` ``axes_ranges`` entry given as ``1D`` array when data is ``>1D``, i.e. as an independent axis. (`#3283 <https://github.com/sunpy/sunpy/pull/3283>`__)
- Fixed a `sunpy.coordinates` bug where a frame using the default observer of Earth could have its observer overwritten during a transformation. (`#3291 <https://github.com/sunpy/sunpy/pull/3291>`__)
- Fixed a bug where the transformation from `~sunpy.coordinates.frames.Helioprojective` to `~sunpy.coordinates.frames.Heliocentric` used the Sun-observer distance from the wrong frame when shifting the origin, and thus might not give the correct answer if the observer was not the same for the two frames. (`#3291 <https://github.com/sunpy/sunpy/pull/3291>`__)
- Fixed a bug with the transformations between `~sunpy.coordinates.frames.Heliocentric` and `~sunpy.coordinates.frames.HeliographicStonyhurst` when the frame observation time was not the same as the observer observation time.  The most common way to encounter this bug was when transforming from `~sunpy.coordinates.frames.Helioprojective` to any non-observer-based frame while also changing the observation time. (`#3291 <https://github.com/sunpy/sunpy/pull/3291>`__)
- VSO client ``fetch`` should not download when ``wait`` keyword argument is specified. (`#3298 <https://github.com/sunpy/sunpy/pull/3298>`__)
- Fixed a bug with `~sunpy.coordinates.wcs_utils.solar_frame_to_wcs_mapping` that assumed that the supplied frame was a SunPy frame. (`#3305 <https://github.com/sunpy/sunpy/pull/3305>`__)
- Fixed bugs with `~sunpy.coordinates.wcs_utils.solar_frame_to_wcs_mapping` if the input frame does not include an observation time or an observer. (`#3305 <https://github.com/sunpy/sunpy/pull/3305>`__)
- `~sunpy.coordinates.utils.GreatArc` now accounts for the start and end points of the arc having different observers. (`#3334 <https://github.com/sunpy/sunpy/pull/3334>`__)
- Fixed situations where 2D coordinates provided to `~sunpy.coordinates.frames.HeliographicStonyhurst` and `~sunpy.coordinates.frames.HeliographicCarrington` were not converted to 3D as intended.  Furthermore, the stored data will always be the post-conversion, 3D version. (`#3351 <https://github.com/sunpy/sunpy/pull/3351>`__)
- Fix off by one error in :func:`sunpy.map.header_helper.make_fitswcs_header` where when using the
  default ``reference_pixel=None`` keyword argument the pixel coordinate of the
  reference pixel was off by +1. (`#3356 <https://github.com/sunpy/sunpy/pull/3356>`__)
- Updated both GOES XRS and LYRA dataretriever clients to use ``sunpy.util.scraper.Scraper``, to make sure that files are actually on the servers being queried. (`#3367 <https://github.com/sunpy/sunpy/pull/3367>`__)
- Fixing the ordering of lon and lat inputs into make_fitswcs_header (`#3371 <https://github.com/sunpy/sunpy/pull/3371>`__)
- Updated the URL for Fermi spacecraft-pointing files to use an HTTPS connection to HEASARC. (`#3381 <https://github.com/sunpy/sunpy/pull/3381>`__)
- Fixed a bug where permission denied errors when downloading files are very verbose by adding an error message in `~sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch`. (`#3417 <https://github.com/sunpy/sunpy/pull/3417>`__)
- Fixed a malformed call to `astropy.time.Time` in a test, which resulted in an incorrect time scale (UTC instead of TT). (`#3418 <https://github.com/sunpy/sunpy/pull/3418>`__)
- Fix incorrect files being included in the tarball, and docs missing from the
  tarball (`#3423 <https://github.com/sunpy/sunpy/pull/3423>`__)
- Fixed a bug where clipping behavior had been enabled by default in the plotting normalizers for ``Map`` objects.  Clipping needs to be disabled to make use of the over/under/masked colors in the colormap. (`#3427 <https://github.com/sunpy/sunpy/pull/3427>`__)
- Fix a bug with observer based frames that prevented a coordinate with an array of obstimes being transformed to other frames. (`#3455 <https://github.com/sunpy/sunpy/pull/3455>`__)
- `sunpy.map.GenericMap` will no longer raise a warning if the position of the
  observer is not known for frames that don't need an observer, i.e. heliographic
  frames. (`#3462 <https://github.com/sunpy/sunpy/pull/3462>`__)
- Apply `os.path.expanduser` to `sunpy.map.map_factory.MapFactory` input
  before passing to `glob.glob` (`#3477 <https://github.com/sunpy/sunpy/pull/3477>`__)
- Fix multiple instances of `sunpy.map` sources assuming the type of FITS Header
  values. (`#3497 <https://github.com/sunpy/sunpy/pull/3497>`__)
- Fixed a bug with `~sunpy.coordinates.NorthOffsetFrame` where non-spherical representations for the north pole produced an error. (`#3517 <https://github.com/sunpy/sunpy/pull/3517>`__)
- Fixed ``map.__repr__`` when the coordinate system information contained in the
  ``CUNIT1/2`` metadata is not set to a known value. (`#3569 <https://github.com/sunpy/sunpy/pull/3569>`__)
- Fixed bugs with some coordinate transformations when ``obstime`` is ``None`` on the destination frame but can be assumed to be the same as the ``obstime`` of the source frame. (`#3576 <https://github.com/sunpy/sunpy/pull/3576>`__)
- Updated `sunpy.map.mapsequence.MapSequence` so that calling ``_derotate()`` raises ``NotImplementedError``.
  Added associated tests. (`#3613 <https://github.com/sunpy/sunpy/pull/3613>`__)
- Fixed pandas plotting registration in `sunpy.timeseries`. (`#3633 <https://github.com/sunpy/sunpy/pull/3633>`__)
- Correctly catch and emit a warning when converting a map metadata to a FITS
  header and it contains a keyword with non-ascii characters. (`#3645 <https://github.com/sunpy/sunpy/pull/3645>`__)


Improved Documentation
----------------------

- Clean up the docstring for `sunpy.physics.differential_rotation.solar_rotate_coordinate` to make the example clearer. (`#2708 <https://github.com/sunpy/sunpy/pull/2708>`__)
- Added new gallery examples and cleaned up various gallery examples. (`#3181 <https://github.com/sunpy/sunpy/pull/3181>`__)
- Cleaned and expanded upon the docstrings for each Fido Client. (`#3220 <https://github.com/sunpy/sunpy/pull/3220>`__)
- Added clarifying hyperlinks to the gallery example ``getting_lasco_observer_location`` to link to `astroquery <https://astroquery.readthedocs.io/en/latest/>`__ docs page. (`#3228 <https://github.com/sunpy/sunpy/pull/3228>`__)
- Added more details to docstrings in `sunpy.coordinates.frames`. (`#3262 <https://github.com/sunpy/sunpy/pull/3262>`__)
- Added a link to package maintainer list in the API Stability page. (`#3281 <https://github.com/sunpy/sunpy/pull/3281>`__)
- Improved the contributing guide by updating commands and highlighting text. (`#3394 <https://github.com/sunpy/sunpy/pull/3394>`__)
- Removing ``.fits`` from the end of path kwargs in `sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch` docs to change output file extension from ``{file}.fits.fits`` to ``{file}.fits``. (`#3399 <https://github.com/sunpy/sunpy/pull/3399>`__)
- A new example gallery section "Using SunPy with Other Packages" has been added,
  which contains a set of new examples using the `reproject
  <https://reproject.readthedocs.io/>`__ with solar data. (`#3405 <https://github.com/sunpy/sunpy/pull/3405>`__)
- Added a table of supported coordinate systems and other miscellaneous improvements to the :ref:`coordinates documentation <sunpy-topic-guide-coordinates-index>`. (`#3414 <https://github.com/sunpy/sunpy/pull/3414>`__)
- Clarified the meaning of :attr:`sunpy.map.GenericMap.dsun`. (`#3430 <https://github.com/sunpy/sunpy/pull/3430>`__)
- Fixed the plots with multiple subplots in the ``Map`` user guide to properly use `~astropy.visualization.wcsaxes` and to be appropriately sized. (`#3454 <https://github.com/sunpy/sunpy/pull/3454>`__)
- Fixed various issues with the gallery example of saving/loading coordinates using ``asdf``. (`#3473 <https://github.com/sunpy/sunpy/pull/3473>`__)
- Added ``sunpy.__citation__`` with a BibTex entry for citing sunpy. (`#3478 <https://github.com/sunpy/sunpy/pull/3478>`__)
- Added an example showing how to display two maps and fade between them. (`#3488 <https://github.com/sunpy/sunpy/pull/3488>`__)
- Clarified the meaning of some `~sunpy.map.GenericMap` observer properties. (`#3585 <https://github.com/sunpy/sunpy/pull/3585>`__)
- Added inherited members of `sunpy.map` classes to the docs. (`#3587 <https://github.com/sunpy/sunpy/pull/3587>`__)
- Fixed documentation of ``sunpy.database.Database.search`` by adding ``Returns`` docstring. (`#3593 <https://github.com/sunpy/sunpy/pull/3593>`__)
- Updated the docstring for the parameter ``sortby`` in `~sunpy.map.MapSequence` with the default value, valid value and how to disable sorting. (`#3601 <https://github.com/sunpy/sunpy/pull/3601>`__)
- Updated the tour guide to reflect that the time series is not random data. (`#3603 <https://github.com/sunpy/sunpy/pull/3603>`__)
- Fixes bold type and extra line breaks of remote data manager example. (`#3615 <https://github.com/sunpy/sunpy/pull/3615>`__)


Trivial/Internal Changes
------------------------

- Allow running our sphinx-gallery examples as Jupyter notebooks via Binder (`#3256 <https://github.com/sunpy/sunpy/pull/3256>`__)
- Improve error messages and type checking in
  ``sunpy.visualization.animator.image.ImageAnimatorWCS``. (`#3346 <https://github.com/sunpy/sunpy/pull/3346>`__)
- Copy the library ``distro`` into :file:`sunpy/extern`: replaces the deprecated ``platform/linux_distribution`` (`#3396 <https://github.com/sunpy/sunpy/pull/3396>`__)
- The version of Matplotlib used to generate figure tests has been bumped from
  3.0.3 to 3.1.1. (`#3406 <https://github.com/sunpy/sunpy/pull/3406>`__)
- Corrected spelling of 'plotting' in timeseries method (changed 'ploting' to 'plotting'). (`#3429 <https://github.com/sunpy/sunpy/pull/3429>`__)
- Switched to "importlib_metadata" to get package version to speed up import of SunPy. (`#3449 <https://github.com/sunpy/sunpy/pull/3449>`__)
- Fix tests for `sunpy.data.data_manager` and ensure they are correctly executed with pytest. (`#3550 <https://github.com/sunpy/sunpy/pull/3550>`__)


1.0.0 (2019-06-01)
==================

Backwards Incompatible Changes
------------------------------

- Move the matplotlib animators from ``sunpy.visualisation.imageanimator`` and
  ``sunpy.visualization.mapcubeanimator`` to `sunpy.visualization.animator`. (`#2515 <https://github.com/sunpy/sunpy/pull/2515>`__)
- Make `sunpy.time.parse_time` return `astropy.time.Time` instead of `datetime.datetime`. (`#2611 <https://github.com/sunpy/sunpy/pull/2611>`__)
- The properties and methods of `sunpy.time.TimeRange` returns `astropy.time.Time` and `astropy.time.TimeDelta` instead of `datetime.datetime` and `datetime.timedelta` respectively. (`#2638 <https://github.com/sunpy/sunpy/pull/2638>`__)
- The ``sunpy.instr.goes`` module now accepts and returns
  `sunpy.timeseries.sources.XRSTimeSeries` objects only. (`#2666 <https://github.com/sunpy/sunpy/pull/2666>`__)
- ``obstime`` keyword param of ``sunpy.instr.goes._goes_lx`` takes a non-scalar `astropy.time.Time` object instead of `numpy.ndarray`. The precision of times contained in `sunpy.timeseries` has been increased to 9 from 6. (`#2676 <https://github.com/sunpy/sunpy/pull/2676>`__)
- Removed ``sunpy.net.jsoc.attrs.Time`` because it served the same purpose as `sunpy.net.attrs.Time` after the switch to `astropy.time.Time`. (`#2694 <https://github.com/sunpy/sunpy/pull/2694>`__)
- Remove unused ``**kwargs`` within TimeSeries functions. (`#2717 <https://github.com/sunpy/sunpy/pull/2717>`__)
- Rotation matrices inside map objects were previously stored as numpy matrices, but are now
  stored as numpy arrays, as numpy will eventually remove their matrix datatype. See
  https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html for more information. (`#2719 <https://github.com/sunpy/sunpy/pull/2719>`__)
- The ``sunpy.cm.show_colormaps`` function now accepts the keyword 'search' instead of 'filter'. (`#2731 <https://github.com/sunpy/sunpy/pull/2731>`__)
- The keyword arguments to all client ``.fetch`` methods have been changed to
  support the new parfive downloader and to ensure consistency across all Fido
  clients. (`#2797 <https://github.com/sunpy/sunpy/pull/2797>`__)
- The Helioviewer client has been switched to using the newer Helioviewer API.
  This has meant that we have changed some of the keywords that were passed into client's methods.
  We have enforced that several keywords (observatory,instrument,detector,measurement) need to be defined otherwise the functions cannot return any data. (`#2801 <https://github.com/sunpy/sunpy/pull/2801>`__)
- Maps no longer assume that the pixel units are arcseconds if the units aren't
  explicitly set. In addition to this if critical metadata is missing from when
  creating a map, the map will fail to initialize and will raise an error. (`#2847 <https://github.com/sunpy/sunpy/pull/2847>`__)
- axis_ranges kwarg of ``sunpy.visualization.animator.base.ArrayAnimator``, ``sunpy.visualization.animator.image.ImageAnimator`` and ``sunpy.visualization.animator.line.LineAnimator`` now must be entered as None, [min, max] or pixel edges of each array element. Previously, pixel centers were expected.  This change removes ambiguity in interpretation and ensures the extent of the plot can always be accurately derived. (`#2867 <https://github.com/sunpy/sunpy/pull/2867>`__)
- All keywords have been added (with defaults) to each ``~sunpy.net.helioviewer.HelioviewerClient`` function.
  This means that there will be some changes to the style of the PNG screenshot that is returned.
  Returns for the JPEG 2000 and the other functions should be the same but not guaranteed. (`#2883 <https://github.com/sunpy/sunpy/pull/2883>`__)
- Changed `sunpy.sun.models.interior` and `sunpy.sun.models.evolution` from `pandas.DataFrame` to `astropy.table.QTable` (`#2936 <https://github.com/sunpy/sunpy/pull/2936>`__)
- Minimum numpy version is now >=1.14.5 (`#2954 <https://github.com/sunpy/sunpy/pull/2954>`__)
- Removed ``sunpy.time.julian_day``, ``sunpy.time.julian_centuries``, ``sunpy.time.day_of_year``, ``sunpy.time.break_time``, ``sunpy.time.get_day``. (`#2999 <https://github.com/sunpy/sunpy/pull/2999>`__)
- Updated the solar values in `sunpy.sun.constants` to IAU 2015 values. (`#3001 <https://github.com/sunpy/sunpy/pull/3001>`__)
- Renamed ``eccentricity_sunearth_orbit`` to ``eccentricity_sun_earth_orbit``. (`#3001 <https://github.com/sunpy/sunpy/pull/3001>`__)
- Renamed ``sunpy.image.rescale`` to `sunpy.image.resample`. (`#3044 <https://github.com/sunpy/sunpy/pull/3044>`__)
- Remove the ``basic_plot`` keyword argument from
  `~sunpy.map.GenericMap.peek`. An example has been added to the gallery
  showing how to make a plot like this. (`#3109 <https://github.com/sunpy/sunpy/pull/3109>`__)
- `sunpy.map.GenericMap` will no longer use the key ``solar_b0`` as a value for heliographic latitude. (`#3115 <https://github.com/sunpy/sunpy/pull/3115>`__)
- `sunpy.map.GenericMap` now checks for a complete observer location rather than
  individually defaulting coordinates (lat, lon, distance) to Earth position. If
  any one of the three coordinates is missing from the header the observer will
  be defaulted to Earth and a warning raised. (`#3115 <https://github.com/sunpy/sunpy/pull/3115>`__)
- ``sunpy.sun.sun`` functions have been re-implemented using Astropy for significantly improved accuracy.  Some functions have been removed. (`#3137 <https://github.com/sunpy/sunpy/pull/3137>`__)
- All of the functions in ``sunpy.sun.sun`` and all of the Sun-specific functions in `sunpy.coordinates.ephemeris` have been moved to the new module `sunpy.coordinates.sun`. (`#3163 <https://github.com/sunpy/sunpy/pull/3163>`__)


Deprecations and Removals
-------------------------

- The deprecated ``sunpy.lightcurve``, ``sunpy.wcs`` and ``sunpy.spectra`` modules have now
  been removed. (`#2666 <https://github.com/sunpy/sunpy/pull/2666>`__)
- ``sunpy.instr.rhessi.get_obssumm_dbase_file`` ``sunpy.instr.rhessi.get_obssum_filename``, ``sunpy.instr.rhessi.get_obssumm_file`` have been removed. `~sunpy.net.Fido` should be used to download these files. (`#2808 <https://github.com/sunpy/sunpy/pull/2808>`__)
- Removed ``heliographic_solar_center`` in favour of ``sunpy.coordinates.get_sun_L0`` and ``sunpy.coordinates.get_sun_B0`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``GenericClient.query`` in favour of `sunpy.net.dataretriever.GenericClient.search` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunearth_distance`` in favour of ``get_sunearth_distance`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``remove_lytaf_events_from_lightcurve`` in favour of ``sunpy.instr.lyra.remove_lytaf_events_from_timeseries`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunpy.cm.get_cmap`` in favour of ``plt.get_cmap`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``database.query`` in favour of ``sunpy.database.Database.search`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunpy.net.vso.InteractiveVSOClient`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``MapCube`` in favour of `~sunpy.map.MapSequence` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``solar_north`` in favour of ``get_sun_P`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``database.download`` in favour of ``sunpy.database.Database.fetch`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunpy.map.GenericMap.pixel_to_data`` in favour of `sunpy.map.GenericMap.pixel_to_world` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``GenericClient.get`` in favour of `sunpy.net.dataretriever.GenericClient.fetch`. This changes applies to the other clients as well. (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``Map.xrange`` and ``Map.yrange`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``sunpy.net.attrs.Wave`` in favour of ``sunpy.net.vso.attrs.Wavelength`` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- Removed ``JSOCClient.check_request`` in favour of `drms.ExportRequest.status` (`#2830 <https://github.com/sunpy/sunpy/pull/2830>`__)
- ``sunpy.net.vso.VSOClient.query_legacy`` and ``sunpy.net.vso.VSOClient.latest`` have been deprecated as we strongly recommend people use `sunpy.net.Fido` for all queries. (`#2866 <https://github.com/sunpy/sunpy/pull/2866>`__)
- The deprecated ``sunpy.physics.transforms`` module has been removed, it is
  replaced by ``sunpy.physics.solar_rotation`` and
  `sunpy.physics.differential_rotation`. (`#2994 <https://github.com/sunpy/sunpy/pull/2994>`__)
- Removed ``sunpy.sun.sun.solar_cycle_number`` because it was fundamentally flawed (`#3150 <https://github.com/sunpy/sunpy/pull/3150>`__)


Features
--------

- Change arguments to ``sunpy.test`` from ``offline=`` and ``online=`` to ``online`` and ``online_only``. This matches the behavior of the figure keyword arguments and comes as a part of a move to using a modified version of the Astropy test runner. (`#1983 <https://github.com/sunpy/sunpy/pull/1983>`__)
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
- The download manager for `~sunpy.net.Fido` has been replaced with
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
- Provided access to the Helioviewer header information using ``~sunpy.net.helioviewer.HelioviewerClient.get_jp2_header`` function. (`#2904 <https://github.com/sunpy/sunpy/pull/2904>`__)
- Add a new WSDL URL and port to support SunPy use of VSO instance at SDAC. (`#2912 <https://github.com/sunpy/sunpy/pull/2912>`__)
- Add support for COSMO K-Coronograph (KCOR) FITS data. (`#2916 <https://github.com/sunpy/sunpy/pull/2916>`__)
- Add logger messaging system based on `~astropy.logger.AstropyLogger`, cleaned up all warnings, removed all print statements. (`#2980 <https://github.com/sunpy/sunpy/pull/2980>`__)
- The function ``sunpy.image.coalignment.get_correlation_shifts`` now issues an error when the number of dimensions
  are not correct instead of a warning and returning None. (`#2980 <https://github.com/sunpy/sunpy/pull/2980>`__)
- The default location of the sunpy sample data has changed to be in the platform
  specific data directory as provided by `appdirs <https://github.com/ActiveState/appdirs>`__. (`#2993 <https://github.com/sunpy/sunpy/pull/2993>`__)
- Add timeseries support for EVE/ESP level 1 data in `sunpy.timeseries.sources` (`#3032 <https://github.com/sunpy/sunpy/pull/3032>`__)
- The default style for Map plots have changed to reflect the changes in Astropy
  3.2. (`#3054 <https://github.com/sunpy/sunpy/pull/3054>`__)
- `sunpy.coordinates.ephemeris.get_body_heliographic_stonyhurst` can now account for light travel time when computing the (apparent) body position, as long as the observer location is provided. (`#3055 <https://github.com/sunpy/sunpy/pull/3055>`__)
- Added a helper function (:func:`sunpy.map.header_helper.make_fitswcs_header`) that allows users to create a meta header for custom created `sunpy.map.GenericMap`. (`#3083 <https://github.com/sunpy/sunpy/pull/3083>`__)
- Map plotting now accepts the optional keyword ``clip_interval`` for specifying a percentile interval for clipping.  For example, if the interval (5%, 99%) is specified, the bounds of the z axis are chosen such that the lowest 5% of pixels and the highest 1% of pixels are excluded. (`#3100 <https://github.com/sunpy/sunpy/pull/3100>`__)
- The new function `~sunpy.coordinates.get_horizons_coord` enables querying JPL HORIZONS for the locations of a wide range of solar-system bodies, including spacecraft. (`#3113 <https://github.com/sunpy/sunpy/pull/3113>`__)


Bug Fixes
---------

- Fix the bug that prevented VSO queries for HMI data from downloading file
  without specifying ``a.Physobs``. (`#2621 <https://github.com/sunpy/sunpy/pull/2621>`__)
- Fix ``sunpy.map.mapcube.MapCube.plot``. The code had not been updated to support the changes to the wcsaxes helper functions. (`#2627 <https://github.com/sunpy/sunpy/pull/2627>`__)
- Replace all use of the deprecated ``sunpy.cm.get_cmap`` with ``matplotlib.cm.get_cmap`` to prevent deprecation warnings being raised. (`#2635 <https://github.com/sunpy/sunpy/pull/2635>`__)
- Fix generation of the coordinate transformation graph with Astropy 3.1.dev (`#2636 <https://github.com/sunpy/sunpy/pull/2636>`__)
- Prevent helioviewer from erroring when downloading file to a directory that
  does not exist. It will now create the directory when required. (`#2642 <https://github.com/sunpy/sunpy/pull/2642>`__)
- Fix transformations into/out of Heliographic Stonyhurst frame when
  the coordinate representation is Cartesian. (`#2646 <https://github.com/sunpy/sunpy/pull/2646>`__)
- Running the figure tests with ``setup.py test`` now saves the figures and the hashes to the same directory as setup.py. (`#2658 <https://github.com/sunpy/sunpy/pull/2658>`__)
- ``sunpy.instr.fermi.met_to_utc`` now returns the correct utc time which takes into account the leap seconds that have passed. (`#2679 <https://github.com/sunpy/sunpy/pull/2679>`__)
- Support passing Python file objects to ``sunpy.io.fits.write``. (`#2688 <https://github.com/sunpy/sunpy/pull/2688>`__)
- Added DRMS to setup.py so sunpy[all] installs it as a dependency. (`#2693 <https://github.com/sunpy/sunpy/pull/2693>`__)
- Fix eve 0cs timeseries separator regex to support Python 3.7 (`#2697 <https://github.com/sunpy/sunpy/pull/2697>`__)
- Fix the bug which crashes `~sunpy.map.sources.LASCOMap` for when 'date-obs' is reformatted again from a self applied function. (`#2700 <https://github.com/sunpy/sunpy/pull/2700>`__)
- Change all instances of quantity_allclose to `astropy.units.allclose` this prevents pytest being needed to import `sunpy.coordinates` on Astropy 3 (`#2701 <https://github.com/sunpy/sunpy/pull/2701>`__)
- Fix RHESSI obssum file downloading to include the final day in the time range. (`#2714 <https://github.com/sunpy/sunpy/pull/2714>`__)
- Raise an error when transforming between HPC and HCC frames if the observer is not the same. (`#2725 <https://github.com/sunpy/sunpy/pull/2725>`__)
- Replaces the existing LASCO C2 and C3 color maps with new ones that perform better with JP2 and Level 0.5, 1 data. (`#2731 <https://github.com/sunpy/sunpy/pull/2731>`__)
- Do not attempt to save a FITS header comment for a keyword which is not in the header. This prevents an error on saving some maps after the metadata had been modified but not the comments. (`#2748 <https://github.com/sunpy/sunpy/pull/2748>`__)
- Add support for `~sunpy.map.sources.HMIMap` objects as input to ``sunpy.instr.aia.aiaprep``. (`#2749 <https://github.com/sunpy/sunpy/pull/2749>`__)
- User can convert between HPC and HCC coordinates with different observers. This is implemented by automatically transforming the coordinate into HGS and then changing observer, and then transforming back to HCC. (`#2754 <https://github.com/sunpy/sunpy/pull/2754>`__)
- Changed default file type for Helioviewer to prevent decode errors. (`#2771 <https://github.com/sunpy/sunpy/pull/2771>`__)
- Increase figure size to avoid cutting off longer colormap names in ``sunpy.cm.show_colormaps``. (`#2824 <https://github.com/sunpy/sunpy/pull/2824>`__)
- The sample data directory will no longer be created until files are downloaded
  to it. (`#2836 <https://github.com/sunpy/sunpy/pull/2836>`__)
- Timeseries and lightcurve will now respect updated config values for download directory. (`#2844 <https://github.com/sunpy/sunpy/pull/2844>`__)
- Always use _default_wrap_angle rather than hard coding a wrap angle in the init
  of a sunpy coordinate frame (`#2853 <https://github.com/sunpy/sunpy/pull/2853>`__)
- Ensure imageanimators only slice arrays with integers (`#2856 <https://github.com/sunpy/sunpy/pull/2856>`__)
- Fixed ``sunpy.io.fits.write`` to handle the keyword ``COMMENT`` correctly. (`#2880 <https://github.com/sunpy/sunpy/pull/2880>`__)
- If Carrington longitude ("crln_obs") is found in the FITS header, `~sunpy.map.Map` converts this to the correct Heliographic longitude. (`#2946 <https://github.com/sunpy/sunpy/pull/2946>`__)
- ``sunpy.net.helio.hec.HECClient.time_query`` now resolves the correct input time format. (`#2969 <https://github.com/sunpy/sunpy/pull/2969>`__)
- Fixes the calculation of the solar rotation of coordinates and the differential rotation of `sunpy.map.GenericMap`. (`#2972 <https://github.com/sunpy/sunpy/pull/2972>`__)
- Added back the FERMI GBM client to `sunpy.net.dataretriever`. (`#2983 <https://github.com/sunpy/sunpy/pull/2983>`__)
- Fix bug in `sunpy.net.hek` which raised and error if a search returned zero results, now returns an empty `sunpy.net.hek.HEKTable`. (`#3046 <https://github.com/sunpy/sunpy/pull/3046>`__)
- `~sunpy.map.sources.AIAMap` now uses the provided HAE coordinates instead of the provided HGS coordinates to determine the observer location. (`#3056 <https://github.com/sunpy/sunpy/pull/3056>`__)
- Correctly zero pad milliseconds in the ``sunpy.util.scraper.Scraper`` formatting to prevent errors when the millisecond value was less than 100. (`#3063 <https://github.com/sunpy/sunpy/pull/3063>`__)
- Fix ``sunpy.util.scraper.Scraper`` failing if a directory is not found on a remote server. (`#3063 <https://github.com/sunpy/sunpy/pull/3063>`__)
- Correctly extract observer location from MDI and EIT data (`#3067 <https://github.com/sunpy/sunpy/pull/3067>`__)
- Fix HGS <> HCRS test due to Ecliptic frame changes in astropy 3.2 (`#3075 <https://github.com/sunpy/sunpy/pull/3075>`__)
- Fixes bug when creating a timeseries from a URL and bug when creating a TimeSeries from  older GOES/XRS fits files. (`#3081 <https://github.com/sunpy/sunpy/pull/3081>`__)
- Added `~sunpy.map.sources.EUVIMap.rsun_obs`. It returns a quantity in arcsec consistent with other `sunpy.map.GenericMap` and overwrites mapbase's assumption of a photospheric limb as seen from Earth. (`#3099 <https://github.com/sunpy/sunpy/pull/3099>`__)
- Fixed bugs related to using `~sunpy.map.GenericMap.plot` and `~sunpy.map.GenericMap.peek` with the ``inline`` Matplotlib backend in Jupyter notebook. (`#3103 <https://github.com/sunpy/sunpy/pull/3103>`__)
- Make a correction to `sunpy.coordinates.wcs_utils.solar_wcs_frame_mapping` so
  that `astropy.wcs.WCS` objects are correctly converted to
  `sunpy.coordinates.frames` objects irrespective of the ordering of the axes. (`#3116 <https://github.com/sunpy/sunpy/pull/3116>`__)
- The `~sunpy.physics.differential_rotation.solar_rotate_coordinate` function returns a coordinate that accounts for the location of the new observer. (`#3123 <https://github.com/sunpy/sunpy/pull/3123>`__)
- Add support for rotation parameters to :func:`sunpy.map.header_helper.make_fitswcs_header`. (`#3139 <https://github.com/sunpy/sunpy/pull/3139>`__)
- Improve the implementation of `~sunpy.physics.differential_rotation.differential_rotate` the image warping when transforming Maps for differential rotation and change in observer position. (`#3149 <https://github.com/sunpy/sunpy/pull/3149>`__)
- Fix a bug where new helioviewer sources potentially cause ``~sunpy.net.helioviewer.HelioviewerClient.data_sources`` to error. (`#3162 <https://github.com/sunpy/sunpy/pull/3162>`__)


Improved Documentation
----------------------

- Organise the gallery into sections based on example type and tidy up a little. (`#2624 <https://github.com/sunpy/sunpy/pull/2624>`__)
- Added gallery example showing the conversion of Helioprojective Coordinates to Altitude/Azimuth Coordinates to and back. (`#2656 <https://github.com/sunpy/sunpy/pull/2656>`__)
- Add contribution guidelines for the sunpy example gallery. (`#2682 <https://github.com/sunpy/sunpy/pull/2682>`__)
- Added a gallery example for "Downloading and plotting a HMI image" and "Creating a Composite map". (`#2746 <https://github.com/sunpy/sunpy/pull/2746>`__)
- Added an example for ``sunpy.visualization.animator.ImageAnimatorWCS``. (`#2752 <https://github.com/sunpy/sunpy/pull/2752>`__)
- Minor changes to the developer guide regarding sprint labels. (`#2765 <https://github.com/sunpy/sunpy/pull/2765>`__)
- Copyedited and corrected the solar cycles example. (`#2770 <https://github.com/sunpy/sunpy/pull/2770>`__)
- Changed "online" mark to "remote_data" and made formatting of marks consistent. (`#2799 <https://github.com/sunpy/sunpy/pull/2799>`__)
- Add a missing plot to the end of the units and coordinates guide. (`#2813 <https://github.com/sunpy/sunpy/pull/2813>`__)
- Added gallery example showing how to access the SunPy colormaps (`#2865 <https://github.com/sunpy/sunpy/pull/2865>`__)
- Added gallery example showing how to access the SunPy solar physics constants. (`#2882 <https://github.com/sunpy/sunpy/pull/2882>`__)
- Major clean up of the developer documentation. (`#2951 <https://github.com/sunpy/sunpy/pull/2951>`__)
- Overhaul of the install instructions for the guide section of our documentation. (`#3147 <https://github.com/sunpy/sunpy/pull/3147>`__)


Trivial/Internal Changes
------------------------

- `~sunpy.time.parse_time` now uses `~functools.singledispatch` underneath. (`#2408 <https://github.com/sunpy/sunpy/pull/2408>`__)
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
- Fix misspelling of private attribute ``_default_heliographic_latitude`` in map. (`#2730 <https://github.com/sunpy/sunpy/pull/2730>`__)
- Miscellaneous fixes to developer docs about building sunpy's documentation. (`#2825 <https://github.com/sunpy/sunpy/pull/2825>`__)
- Changed ``sunpy.instr.aia.aiaprep`` to update BITPIX keyword to reflect the float64 dtype. (`#2831 <https://github.com/sunpy/sunpy/pull/2831>`__)
- Remove warning from ``GenericMap.submap`` when using pixel ``Quantities`` as input. (`#2833 <https://github.com/sunpy/sunpy/pull/2833>`__)
- Remove the usage of six and all ``__future__`` imports (`#2837 <https://github.com/sunpy/sunpy/pull/2837>`__)
- Fix SunPy Coordinate tests with Astropy 3.1 (`#2838 <https://github.com/sunpy/sunpy/pull/2838>`__)
- Stores entries from directories into database sorted by name. It adds mocks to the database user guide examples. (`#2873 <https://github.com/sunpy/sunpy/pull/2873>`__)
- Fix all DeprecationWarning: invalid escape sequence. (`#2885 <https://github.com/sunpy/sunpy/pull/2885>`__)
- Used `unittest.mock` for creating offline tests for simulating online tests for :file:`test_noaa.py` (`#2900 <https://github.com/sunpy/sunpy/pull/2900>`__)
- Fix support for pip 19 and isolated builds (`#2915 <https://github.com/sunpy/sunpy/pull/2915>`__)
- Moved to using `AppDirs <https://github.com/ActiveState/appdirs>`__ as the place to host our configuration file. (`#2922 <https://github.com/sunpy/sunpy/pull/2922>`__)
- Users can now use fewer keywords in our ``~sunpy.net.helioviewer.HelioviewerClient`` to access the available sources. Either by ``observatory`` and ``measurement`` or ``instrument`` and ``measurement`` as this much information is enough to get the source ID for most of the cases. (`#2926 <https://github.com/sunpy/sunpy/pull/2926>`__)
- Remove the pytest dependency on the ``GenericMap`` asdf tag. (`#2943 <https://github.com/sunpy/sunpy/pull/2943>`__)
- Fix initialization of `~sunpy.net.vso.VSOClient` when no WSDL link is found. (`#2981 <https://github.com/sunpy/sunpy/pull/2981>`__)


0.9.0
=====

New Features
------------

- Added TimeUTime class to support utime. [#2409]
- Example for fine-grained use of ticks and grids [#2435]
- Maintiners Workflow Guide [#2411]
- Decorator to append and/or prepend doc strings [#2386]
- Adding ``python setup.py test --figure-only`` [#2557]
- Fido.fetch now accepts pathlib.Path objects for path attribute.[#2559]
- The `~sunpy.coordinates.HeliographicStonyhurst` coordinate system can now be specified
  using a cartesian system, which is sometimes known as the
  "Heliocentric Earth equatorial" (HEEQ) coordinate system. [#2437]

API Changes
-----------

- ``sunpy.coordinates.representation`` has been removed. Longitude wrapping is now done in the constructor of the frames. [#2431]
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
- Added docstrings to functions in download.py [#2415]
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
- removed ``wavelnth`` keyword in meta desc of Maps to avoid using non standard FITS keyword like ``nan`` [#2456]
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
  ``1*u.one``. This fixes a plotting bug with ``WCSAxes`` in Astropy 3.0 [#2465]
- removed ``wavelnth`` keyword in meta desc of Maps to avoid using non standard FITS keyword like ``nan`` [#2427]
- Change the default units for HPC distance from ``u.km`` to `None`. [#2465]

0.8.3
=====

Bug Fixes
---------

- `~sunpy.net.dataretriever.XRSClient` now reports time ranges of files correctly. [#2364]
- Make parse_time work with datetime64s and pandas series [#2370]
- CompositeMap axes scaling now uses map spatial units [#2310]
- Moved license file to root of repository and updated README file [#2326]
- Fix docstring formatting for net.vso.attrs [#2309]]
- Fix coloring of ticks under matplotlib 2.0 default style [#2320]
- Always index arrays with tuples in ``ImageAnimator`` [#2320]
- Added links to possible attrs for FIDO in guide [#2317] [#2289]
- Updated GitHub Readme [#2281] [#2283]
- Fix matplotlib / pandas 0.21 bug in examples [#2336]
- Fixes the off limb enhancement example [#2329]
- Changes to masking hot pixels and picking bright pixels examples [#2325] [#2319]
- Travis CI fix for numpy-dev build [#2340]
- Updated masking brightest pixel example [#2338]
- Changed TRAVIS cronjobs [#2338]
- Support array values for ``obstime`` for coordinates and transformations [#2342] [#2346]
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
   uniform API that is ``search`` and ``fetch``. The older functions are
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

-  Fix test failure (mapbase) with 1.7.4
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
-  HELIO URL updated, queries should now work as expected.
-  All tabs removed from the code base.
-  All tests now use tempfile rather than creating files in the current
   directory.
-  Documentation builds under newer sphinx versions.
-  ANA and JP2 tests are skipped if dependencies are missing.
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
-  Added t-'now' to parse\_time to provide utcnow datetime.
-  Fixed time dependent functions (.sun) to default to t-'now'
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
-  fix file paths to use os.path.join for platform independence.

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
-  Fix apparent Right Ascension calculations.
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
-  GOESLightCurve now fails politely if no data is available.

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

-  Removal of Optional PIL dependency
-  Parse\_time now looks through nested lists/tuples
-  Draw\_limb and draw\_grid are now implemented on MapCube and
   CompositeMap
-  Calculations for differential rotation added
-  mapcube.plot() now runs a mpl animation with optional controls
-  A basic Region of Interest framework now exists under sunpy.roi
-  STEREO COR colour maps have been ported from solarsoft.
-  sunpy.time.timerange has a split() method that divides up a time
   range into n equal parts.
-  Added download progress bar
-  pyfits is deprecated in favor of Astropy

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
groundwork for GenericMap inheriting from astropy.NDData has been done,
there is now a NDDataStandin class to provide basic functionality.

io: \* top level file\_tools improved to be more flexible and support
multiple HDUs \* all functions in sunpy.io now assume multiple HDUs,
even JP2 ones. \* there is now a way to override the automatic filetype
detection \* Automatic fits file detection improved \* extract\_waveunit
added to io.fits for detection of common ways of storing wavelength unit
in fits files.

-  A major re-work of all internal imports has resulted in a much cleaner
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
