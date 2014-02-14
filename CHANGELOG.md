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
Bug Fix: Fix a regression in CompositeMap that made contor plots fail.
Bug Fix: Allow Map() to accept dict as metadata.
Bug Fix: Pass arguments from Map() to io.read_file.

0.3.0
=====
Major Changes:

    Removal of Optional PIL dependancy
    Parse_time now looks through nested lists/tuples
    Draw_limb and draw_grid are now implemented on MapCube and CompositeMap
    Caculations for differential roation added
    mapcube.plot() now runs a mpl animation with optional controls
    A basic Region of Interest framework now exists under sunpy.roi
    STEREO COR colour maps have been ported from solarsoft.
    sunpy.time.timerange has a split() method that divides up a time range into n equal parts.
    Added download progress bar
    pyfits is depricated in favor of Astropy
    
spectra:
    Plotting has been refactorted to use a consistent interface
    spectra now no-longer inherits from numpy.ndarray instead has a .data attribute.

Map:
    map now no-longer inherits from numpy.ndarray instead has a .data attribute.
    make_map is deprecated in favor of Map which is a new factory class
    sunpy.map.Map is now sunpy.map.GenericMap
    mymap.header is now mymap.meta
    attributes of the map class are now read only, changes have to be made through map.meta
    new MapMeta class to replace MapHeader, MapMeta is not returned by sunpy.io
    The groundwork for GenericMap inherting from astropy.NDData has been done,
        there is now a NDDataStandin class to provide basic functionality.
    
io:  
    top level file_tools improved to be more flexible and support multiple HDUs
    all functions in sunpy.io now assume mutliple HDUs, even JP2 ones.
    there is now a way to override the automatic filetype detection
    Automatic fits file detection improved
    extract_waveunit added to io.fits for detection of common ways of storing
        wavelength unit in fits files.
      

Bug fixes or under the hood changes:

    A major re-work of all interal imports has resulted in a much cleaner namespace, i.e. sunpy.util.util is no longer used to import util.
    Some SOHO and STEREO files were not reading properly due to a date_obs parameter.
    Sunpy will now read JP2 files without a comment parameter.
    Memory leak in Crotate patched
    Callisto: Max gap between files removed

0.2.0
=====
Below are the main features that have been added for this release:

Completely re-written plotting routines for most of the core datatypes.
JPEG 2000 support as an input file type.
Improved documentation for much of the code base, including re-written installation instructions.
New lightcurve object
    LYRA support
    GOES/XRS support
    SDO/EVE support
New Spectrum and Spectrogram object (in development)
    Spectrogram plotting routines
    Callisto spectrum type and support
    STEREO/SWAVES support
Map Object
    Added support for LASCO, Yohkoh/XRT maps
    A new CompositeMap object for overlaying maps
    Resample method
    Superpixel method
    The addition of the rotate() method for 2D maps.
