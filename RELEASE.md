The SunPy project is happy to announce the release of SunPy 0.5.0.
This release consists of x commits from y people, including the ability to
co-align map cubes via template matching in scikit image, massive improvements
to Map.rotate() including an implementation of an aiaprep calibration routine for
SDO/AIA data and the ability to calculate GOES temperature and emission
measure from GOES fluxes.

Special mentions to Daniel Ryan and Andrew Leonard who have both contributed
their first major features to SunPy for this release. They contributed the
GOES temperature and aiaprep code respectively.

New Features:

    * map.rotate() improvements and the additon of a aiaprep routine.
    * GOES temperature and emission measure calculation.
    * Added functions that implement image coalignment with support for MapCubes.
    * MapCube._maps changed to MapCube.maps.
    * Added Nobeyama Radioheliograph data support to Lightcurve object.
    * Added support for NOAA solar cycle prediction in lightcurves.
    * Improved line quality and performances issues with map.draw_grid().
    * Most tests should pass on windows and on installed versions of SunPy.
    * Added a window/split method to time range.
    * Updates to spectrogram documentation.
    * Added method Database.add_from_hek_query_result to HEK database.
    * Added method Database.download_from_vso_query_result.
    * GOES Lightcurve now makes use of a new source of GOES data, provides metadata, and data back to 1981.
    * Fix algorithm in sunpy.sun.equation_of_center.
    * Added contains functionality to TimeRange module
    * Added t='now' to parse_time to privide utcnow datetime.
    * Fixed time dependant functions (.sun) to default to t='now'
    * Fixed solar_semidiameter_angular_size
    * Removed sqlalchemy as a requirement for SunPy
    * Some basic tests for GenericLightCurve on types of expected input.
    * Added Docstrings to LightCurve methods.
    * Added tests for classes in sunpy.map.sources.
    * Cleaned up the sunpy namespace, removed .units, /ssw and .sphinx. Also moved .coords .physics.transforms.

The people who have contributed to this release are:

    Stuart Mumford
    Simon Liedtke
    Steven Christe
    Jack Ireland
    Andrew Inglis
    Nabil Freij
    Samuel Bennett
    David Perez-Suarez
    Pritish Chakraborty
    Albert Y. Shih
    John Evans
    Michael Malocha
    Florian Mayer
    Russell Hewett
    Jose Iván Campos
    Keith Hughitt
    Tiago Pereira
