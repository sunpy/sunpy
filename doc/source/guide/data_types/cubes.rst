=====
Cubes
=====

SunPy Cubes are 3- or 4-dimensional arrays of continuous data with an associated coordinate system.
There are four possible cube types:: time and two solar dimensions (i.e. a time-series of maps); energy and two solar dimensions (i.e. a series of maps ordered by wavelength); time, energy and a single solar dimension (i.e. a timeseries cube); or a time, energy and two solar dimensions hypercube.
It is important to note that regardless of the data source used, the data held in the cube is ordered by the present dimension type:: the highest priority is assigned to time, then energy, and solar coordinates at the bottom. This means that if, for example, a cube is given solar-x, energy and solar-y as its dimensions, it will rearrange the axes and data to make it an energy-x-y cube. When more than one spatial dimansion is present, the oriiginal order is preserved.
This document is meant merely as an introduction; for the full documentation consult the API Reference.

------------
Data Support
------------
The cube datatype currently only supports EIS. Future datasources will be added.

1. Creating en EIS Cube
-----------------------
First of all, we need to initialize the current interactive python shell::

    import sunpy
    import sunpy.cube
    from sunpy.cube.sources.eis import  EISSpectralCube as eis
    my_cubes = eis.read("/directory/file.fits")

This is assuming your FITS file's address is /directory/file.fits. EIS convention is to include several cubes in a single file, ordered by their principal wavelength. Therefore, my_cubes is a variable containing a python dictionary of cubes ordered by their principal wavelength.
Let's see a typical example::

    In [5]: my_cubes.keys()
    Out[5]: ['FE XV 284.160',
             'SI VII 275.350',
             'MG VI 270.400',
             'FE XII 195.120',
             'FE X 184.540',
             'CA XVII 192.820',
             'HE II 256.320',
             'FE XIII 202.040',
             'FE XII 186.880',
             'FE XIII 203.830']
    In [6]: my_cube = my_cubes['FE X 184.540']

The variable my_cube is now an EIS Spectral Cube.

2. Creating a Generic Cube
--------------------------
If you have your data and an associated WCS object, you can call the Cube constructor by doing::

    my_generic_cube = Cube(data=my_data, wcs=my_wcs)

The internal rearrangement of data and coordinate system is done independently of the variables you put in, so they are not changed.
Note that, due to the way in which WCS objects are handled, wcs cannot be an astropy wcs object; it must be a SunPy WCS one.

3. Converting to other SunPy Datatypes
--------------------------------------
Cubes can be automatically converted to all other SunPy datatypes either pixel-wise or by using world coordinates.
There are two ways of achieving this:: Using the dedicated methods or by slicing the cube as if it were a numpy array. The dedicated methods, which are all called some variation of slice_to_map, slice_to_lightcurve, etc., are more flexible and powerful but at the same time more verbose. Consult the documentation for specifics.
Continuing from the example above, my_cube is a wavelength-solar y-solar x cube. This means that a slice along the energy axis yields an x-y slice - that's a map! Here's an example of one::

    from astropy import units as u
    my_map = my_cube[184.552 * u.Angstrom]
    my_map.peek()

Numpy-style slicing supports almost the same syntax that numpy does, except using None as an index to create another axis. Ranges and increments are fully supported.
Note that you only have to provide the unit once per axis. That is::

    my_cube[184.540 * u.Angstrom:184.550:0.05]

is just as valid as::

    my_cube[184.540 * u.Angstrom:184.550 * u.Angstrom:0.05 * u.Angstrom]

However, you should keep in mind that units cannot be mixed, even between units that represent the same dimension.
The following table shows how slices are interpreted.

+------------+-----------------------------------+
|  Datatype  |          Resulting Slice          |
+============+===================================+
|    Map     |         Celestial x and y         |
+------------+-----------------------------------+
| Lightcurve |- Time (wavelength present in cube)|
|            |                                   |
|            |- Wavelength and a range of time   |
+------------+-----------------------------------+
| Spectrogram|         Time and Energy           |
+------------+-----------------------------------+
|  Spectrum  |   Spatial or energy and spatial   |
+------------+-----------------------------------+
|    Cube    |          Any 3D result            |
+------------+-----------------------------------+
| Numpy array|        Any other slicing          |
+------------+-----------------------------------+

4. Animating a cube
-------------------
You can view a cube like this::

    my_cube.animate()

This creates a new window with a slider that you can adjust to travel between the different 2D slices.