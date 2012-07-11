-------------
Using the HEK
-------------

The Heliophysics Event Knowledgebase (HEK) is a repository of feature and
event information concerning the Sun.  Entries are generated both by automated
algorithms and human observers.  SunPy accesses this information through the
'hek' module, which was developed through support from the European Space
Agency  Summer of Code in Space (ESA-SOCIS) 2011.

1. Creating maps
----------------
SunPy contains a number of example fits files. To make things easy,
SunPy includes several example files which are used throughout the docs. These
files have names like `sunpy.AIA_171_IMAGE` and `sunpy.RHESSI_IMAGE`.
To create the sample AIA map type the following into your interactive Python shell::

    import sunpy
    my_map = sunpy.make_map(sunpy.AIA_171_IMAGE)

The variable my_map is a SunPy map. To create a map from a local fits file try
something like the following ::

    my_map = sunpy.make_map('/mydirectory/mymap.fits')

SunPy should automatically detect the type of file (e.g. fits), what instrument it is 
associated with (e.g. AIA, EIT, LASCO) and will automatically look in the appropriate places for the fits
keywords it needs to interpret the coordinate system. If the type of fits file 
is not recognized then SunPy will try some default fits keywords but results
may vary. The list of files which are currently supported by SunPy can be found by looking at the 
documentation for make_map(). SunPy can also create maps from the jpg2000 files from
`helioviewer.org <http://helioviewer.org/>`.