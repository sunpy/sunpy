Latest
------
* Fixed `aiaprep` to return properly sized map.

* Deprecation warnings fixed when using image coalignment.
* Sunpy is now Python 3.x compatible (3.4 and 3.5).
* Added a unit check and warnings for map metadata.
* Added IRIS SJI color maps.
* Updated `show_colormaps()` with new string filter to show a subset of color maps.
* Fixed MapCube animations by working around a bug in Astropy's ImageNormalize
* Remove ``vso.QueryResponse.num_records()`` in favour of `len(qr)`
* Add a `draw_rectangle` helper to `GenericMap` which can plot rectangles in the
  native coordinate system of the map.
* Added the ability to shift maps to correct for incorrect map location, for example.
* Bug fix for RHESSI summary light curve values.
* Mapcube solar derotation and coalignment now pass keywords to the routine used to
  shift the images, scipy.ndimage.interpolation.shift.
* Add automatic registration of ``GenericMap`` subclasses with the factory as
  long as they define an ``is_datasource_for`` method.
* Added functions ``flareclass_to_flux`` and ``flux_to_flareclass`` which convert
  between GOES flux to GOES class numbers (e.g. X12, M3.4).
* Removed old ``sunpy.util.goes_flare_class()``
* Bug fix for RHESSI summary light curve values.
* The ``MapCube.as_array`` function now returns a masked numpy array if at least
  one of the input maps in the MapCube has a mask.
* Map superpixel method now respects maps that have masks.
* Map superpixel method now accepts numpy functions as an argument, or any user-defined
  function.
* Map superpixel method no longer has the restriction that the number of original pixels
  in the x (or y) side of the superpixel exactly divides the number of original
  pixels in the x (or y) side of the original map data.
* `sunpy.physics.transforms` has been deprecated and the code moved into `sunpy.physics`.
* Add the `sunpy.coordinates` module, this adds the core physical solar coordinates frame within the astropy coordinates framework.
* Added ability of maps to draw contours on top of themselves (`draw_contours`)
* Added concatenate functionality to lightcurve base class.
* Fix Map to allow astropy.io.fits Header objects as valid input for meta arguments.
* Added an examples gallery using `sphinx-gallery`.
* API clean up to constants. Removed constant() function which is now replaced by get().
* Prevent helioviewer tests from checking access to the API endpoint when running tests offline.


* `GenericMap.units` is renamed to `GenericMap.spatial_units` to avoid confusion with `NDData.unit`.
0.6.0
-----

 * Enforced the use of Astropy Quantities through out most of SunPy.
 * Dropped Support for Python 2.6.
 * Remove old style string formatting and other 2.6 compatibility lines.
 * Added vso like querying feature to JSOC Client.
 * Refactor the JSOC client so that it follows the .query() .get() interface of VSOClient and UnifedDownloader.
 * Provide `__str__` and `__repr__` methods on vso `QueryResponse` deprecate `.show()`.
 * Downloaded files now keep file extensions rather than replacing all periods with underscores.
 * Update to TimeRange API, removed t1 and t0, start and end are now read-only attributes.
 * Added ability to download level3 data for lyra Light Curve along with corresponding tests.
 * Added support for gzipped FITS files.
 * Add STEREO HI Map subclass and color maps.
 * Map.rotate() no longer crops any image data.
 * For accuracy, default Map.rotate() transformation is set to bi-quartic.
 * `sunpy.image.transform.affine_transform` now casts integer data to float64 and sets NaN values to 0 for all transformations except scikit-image rotation with order <= 3.
 * CD matrix now updated, if present, when Map pixel size is changed.
 * Removed now-redundant method for rotating IRIS maps since the functionality exists in Map.rotate()
 * Provide `__str__` and `__repr__` methods on vso `QueryResponse` deprecate `.show()`
 * SunPy colormaps are now registered with matplotlib on import of `sunpy.cm`
 * `sunpy.cm.get_cmap` no longer defaults to 'sdoaia94'
 * Added database url config setting to be setup by default as a sqlite database in the sunpy working directory
 * Added a few tests for the sunpy.roi module
 * Added capability for figure-based tests
 * Removed now-redundant method for rotating IRIS maps since the functionality exists in Map.rotate().
 * SunPy colormaps are now registered with matplotlib on import of `sunpy.cm`.
 * `sunpy.cm.get_cmap` no longer defaults to 'sdoaia94'.
 * Added database url config setting to be setup by default as a sqlite database in the sunpy working directory.
 * Added a few tests for the sunpy.roi module.
 * Refactored mapcube co-alignment functionality.
 * Removed sample data from distribution and added ability to download sample files
 * Changed start of GOES 2 operational time range back to 1980-01-04 so data from 1980 can be read into GOESLightCurve object
 * Require JSOC request data calls have an email address attached.
 * Calculation of the solar rotation of a point on the Sun as seen from Earth, and its application to the de-rotation of mapcubes.
 * Downloaded files now keep file extensions rather than replacing all periods with underscores
 * Fixed the downloading of files with duplicate names in sunpy.database
 * Removed sample data from distribution and added ability to download sample files.
 * Added the calculation of the solar rotation of a point on the Sun as seen from Earth, and its application to the de-rotation of mapcubes.
 * Changed default for GOESLightCurve.create() so that it gets the data from the most recent existing GOES fits file.
 * Map plot functionality now uses the mask property if it is present, allowing the plotting of masked map data
 * Map Expects Quantities and returns quantities for most parameters.
 * Map now used Astropy.wcs for world <-> pixel conversions.
 * map.data_to_pixel now has a similar API to map.pixel_to_data.
 * map.shape has been replaced with map.dimensions, which is ordered
   x first.
 * map.rsun_arcseconds is now map.rsun_obs as it returns a quantity.
 * Map properties are now named tuples rather than dictionaries.
 * Improvement for Map plots, standardization and improved color tables,
   better access to plot variables through new plot_settings variable.
 * Huge improvements in Instrument Map doc strings. Now contain instrument
   descriptions as well as reference links for more info.
 * net.jsoc can query data series with time sampling by a Sample attribute implemented in vso.
 * MapCube.plot and MapCube.peek now support a user defined plot_function argument for customising the animation.
 * Added new sample data file, an AIA cutout file.
 * Moved documentation build directory to doc/build

0.5.0
-----

 * Added additional functionality to the GOES module i.e. the ability to calculate GOES temperature and emission measure from GOES fluxes.
 * changed _maps attribute in MapCube to a non-hidden type
 * Added Nobeyama Radioheliograph data support to Lightcurve object.
 * Fixed some tests on map method to support Windows
 * Added a window/split method to time range
 * Updates to spectrogram documentation
 * Added method Database.add_from_hek_query_result to HEK database
 * Added method Database.download_from_vso_query_result
 * GOES Lightcurve now makes use of a new source of GOES data, provides metadata, and data back to 1981.
 * Removed sqlalchemy as a requirement for SunPy
 * Added support for NOAA solar cycle prediction in lightcurves
 * Some basic tests for GenericLightCurve on types of expected input.
 * Fix algorithm in sunpy.sun.equation_of_center
 * Added Docstrings to LightCurve methods.
 * Added tests for classes in sunpy.map.sources. Note that some classes (TRACE, RHESSI) were left out because SunPy is not able to read their FITS files.
 * Added functions that implement image coalignment with support for MapCubes.
 * Cleaned up the sunpy namespace, removed .units, /ssw and .sphinx. Also moved .coords .physics.transforms.
 * Added contains functionality to TimeRange module
 * Added t='now' to parse_time to privide utcnow datetime.
 * Fixed time dependant functions (.sun) to default to t='now'
 * Fixed solar_semidiameter_angular_size
 * Improved line quality and performances issues with map.draw_grid()
 * Remove deprecated make_map command.

0.4.1
-----
Bug Fixes:
    * Fix map.rotate() functionality
    * Change of source for GOES data.
    * Fix EIT test data and sunpy FITS saving
    * Some documentation fixes
    * fix file paths to use os.path.join for platform independance.


0.4.0
-----
Features:

 * **Major** documentation refactor. A far reaching re-write and restructure.
 * Add a SunPy Database to store and search local data.
 * Add beta support for querying the HELIO HEC
 * Add beta HEK to VSO query translation.
 * Add the ability to download the GOES event list.
 * Add support for downloading and querying the LYTAF database.
 * Add support for ANA data.
 * Updated sun.constants to use astropy.constants objects which include units, source,
 and error instide. For more info check out http://docs.astropy.org/en/latest/constants/index.html
 * Add some beta support for IRIS data products
 * Add a new MapCubeAnimator class with interactive widgets which is returned by mapcube.peek().
 * The Glymur library is now used to read JPEG2000 files.
 * GOESLightCurve now supports all satellites.

Bug Fixes:

 * Add support for VSO queries through proxies.
 * Fix apparent Right Ascension calulations.
 * LightCurve meta data member now an OrderedDict Instance

0.3.2
-----
Bug Fixes:

 * Pass draw_limb arguments to patches.Circle
 * Pass graw_grid arguments to pyplot.plot()
 * Fix README code example
 * Fix Documentation links in potting guide
 * Update to new EVE data URL
 * Update LogicalLightcurve example in docs
 * Improved InteractiveVSOClient documentation
 * GOESLightCurve now fails politely if no data is avalible.

Known Bugs:

 * sunpy.util.unit_conversion.to_angstrom does not work if 'nm' is passed in.

0.3.1
-----

* Bug Fix: Fix a regression in CompositeMap that made contor plots fail.
* Bug Fix: Allow Map() to accept dict as metadata.
* Bug Fix: Pass arguments from Map() to io.read_file.

0.3.0
=====
Major Changes:

 * Removal of Optional PIL dependancy
 * Parse_time now looks through nested lists/tuples
 * Draw_limb and draw_grid are now implemented on MapCube and CompositeMap
 * Caculations for differential roation added
 * mapcube.plot() now runs a mpl animation with optional controls
 * A basic Region of Interest framework now exists under sunpy.roi
 * STEREO COR colour maps have been ported from solarsoft.
 * sunpy.time.timerange has a split() method that divides up a time range into n equal parts.
 * Added download progress bar
 * pyfits is depricated in favor of Astropy

spectra:

 * Plotting has been refactorted to use a consistent interface
 * spectra now no-longer inherits from numpy.ndarray instead has a .data attribute.

Map:
 * map now no-longer inherits from numpy.ndarray instead has a .data attribute.
 * make_map is deprecated in favor of Map which is a new factory class
 * sunpy.map.Map is now sunpy.map.GenericMap
 * mymap.header is now mymap.meta
 * attributes of the map class are now read only, changes have to be made through map.meta
 * new MapMeta class to replace MapHeader, MapMeta is not returned by sunpy.io
 * The groundwork for GenericMap inherting from astropy.NDData has been done, there is now a NDDataStandin class to provide basic functionality.

io:
 * top level file_tools improved to be more flexible and support multiple HDUs
 * all functions in sunpy.io now assume mutliple HDUs, even JP2 ones.
 * there is now a way to override the automatic filetype detection
 * Automatic fits file detection improved
 * extract_waveunit added to io.fits for detection of common ways of storing wavelength unit in fits files.


Bug fixes or under the hood changes:

 * A major re-work of all interal imports has resulted in a much cleaner namespace, i.e. sunpy.util.util is no longer used to import util.
 * Some SOHO and STEREO files were not reading properly due to a date_obs parameter.
 * Sunpy will now read JP2 files without a comment parameter.
 * Memory leak in Crotate patched
 * Callisto: Max gap between files removed

0.2.0
=====
Below are the main features that have been added for this release:

* Completely re-written plotting routines for most of the core datatypes.
* JPEG 2000 support as an input file type.
* Improved documentation for much of the code base, including re-written installation instructions.
* New lightcurve object
    * LYRA support
    * GOES/XRS support
    * SDO/EVE support
* New Spectrum and Spectrogram object (in development)
    * Spectrogram plotting routines
    * Callisto spectrum type and support
    * STEREO/SWAVES support
* Map Object
    * Added support for LASCO, Yohkoh/XRT maps
    * A new CompositeMap object for overlaying maps
    * Resample method
    * Superpixel method
    * The addition of the rotate() method for 2D maps.
