.. _how_to_custom_maps:

Create Custom Maps
==================

It is also possible to create Maps using custom data (e.g. from a simulation or an observation from a data source that is not explicitly supported in **sunpy**).
To do this, you need to provide `sunpy.map.Map` with both the data array as well as appropriate meta information.
The meta information informs `sunpy.map.Map` of the correct coordinate information associated with the data array and should be provided to `sunpy.map.Map` in the form of a header as a `dict` or `~sunpy.util.MetaDict`.
See this :ref:`sphx_glr_generated_gallery_map_map_from_numpy_array.py` for a brief demonstration of generating a Map from a data array.

The keys required for the header information follow the `FITS standard <https://fits.gsfc.nasa.gov/fits_dictionary.html>`__.
**sunpy** provides a Map header helper function to assist in creating a header that contains the correct meta information.
This includes a :func:`~sunpy.map.header_helper.meta_keywords` function that will return a `dict` of the meta keywords used when creating a Map.

.. code-block:: python

    >>> from sunpy.map.header_helper import meta_keywords

    >>> meta_keywords() # doctest: +SKIP
    {'cunit1': 'Units of the coordinate increments along naxis1 e.g. arcsec **required',
     'cunit2': 'Units of the coordinate increments along naxis2 e.g. arcsec **required',
     'crval1': 'Coordinate value at reference point on naxis1 **required'
     ...

The utility function :func:`~sunpy.map.header_helper.make_fitswcs_header` will return a header with the appropriate FITS keywords once the Map data array and an `astropy.coordinates.SkyCoord` or `sunpy.coordinates.frames` is provided.
All the metadata keywords that a Map will parse along with their description are listed in the :ref:`Meta Keywords Table` at the end of this page.

The `astropy.coordinates.SkyCoord` is defined by the user and contains information on the reference frame, reference coordinate, and observer location.
This function returns a `sunpy.util.MetaDict`.
The `astropy.coordinates.SkyCoord` or `sunpy.coordinates.frames` must contain an observation time.

The :func:`~sunpy.map.header_helper.make_fitswcs_header` function also takes optional keyword arguments including ``reference_pixel`` and ``scale`` that describe the pixel coordinate at the reference coordinate (defined by the `~astropy.coordinates.SkyCoord`) and the spatial scale of the pixels, respectively.
If neither of these are given their values default to the center of the data array and 1 arcsec, respectively.

Here's an example of creating a header from some generic data and an `astropy.coordinates.SkyCoord`:

.. code-block:: python

    >>> import numpy as np
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u

    >>> from sunpy.coordinates import frames
    >>> from sunpy.map.header_helper import make_fitswcs_header

    >>> data = np.arange(0,100).reshape(10,10)
    >>> coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime = '2013-10-28', observer = 'earth', frame = frames.Helioprojective)
    >>> header = make_fitswcs_header(data, coord)
    >>> for key, value in header.items():
    ...     print(f"{key}: {value}")
    wcsaxes: 2
    crpix1: 5.5
    crpix2: 5.5
    cdelt1: 1.0
    cdelt2: 1.0
    cunit1: arcsec
    cunit2: arcsec
    ctype1: HPLN-TAN
    ctype2: HPLT-TAN
    crval1: 0.0
    crval2: 0.0
    lonpole: 180.0
    latpole: 0.0
    mjdref: 0.0
    date-obs: 2013-10-28T00:00:00.000
    rsun_ref: 695700000.0
    dsun_obs: 148644585949.49
    hgln_obs: 0.0
    hglt_obs: 4.7711570596394
    naxis: 2
    naxis1: 10
    naxis2: 10
    pc1_1: 1.0
    pc1_2: -0.0
    pc2_1: 0.0
    pc2_2: 1.0
    rsun_obs: 965.3829548285768

From this we can see now that the function returned a `sunpy.util.MetaDict` that populated the standard FITS keywords with information provided by the passed `astropy.coordinates.SkyCoord`, and the data array.
Since the ``reference_pixel`` and keywords were not passed in the example above, the values of ``crpix`` and ``cdelt`` were set to the default values.

These keywords can be passed to the function in the form of an `astropy.units.Quantity` with associated units.
Here's another example of passing ``reference_pixel`` and ``scale`` to the function:

.. code-block:: python

    >>> header = make_fitswcs_header(data, coord,
    ...                                        reference_pixel=u.Quantity([5, 5]*u.pixel),
    ...                                        scale=u.Quantity([2, 2] *u.arcsec/u.pixel))
    >>> for key, value in header.items():
    ...     print(f"{key}: {value}")
    wcsaxes: 2
    crpix1: 6.0
    crpix2: 6.0
    cdelt1: 2.0
    cdelt2: 2.0
    cunit1: arcsec
    cunit2: arcsec
    ctype1: HPLN-TAN
    ctype2: HPLT-TAN
    crval1: 0.0
    crval2: 0.0
    lonpole: 180.0
    latpole: 0.0
    mjdref: 0.0
    date-obs: 2013-10-28T00:00:00.000
    rsun_ref: 695700000.0
    dsun_obs: 148644585949.49
    hgln_obs: 0.0
    hglt_obs: 4.7711570596394
    naxis: 2
    naxis1: 10
    naxis2: 10
    pc1_1: 1.0
    pc1_2: -0.0
    pc2_1: 0.0
    pc2_2: 1.0
    rsun_obs: 965.3829548285768

As we can see, a list of WCS and observer meta information is contained within the generated headers, however we may want to include other meta information including the observatory name, the wavelength and waveunit of the observation.
Any of the keywords in the dictionary returned by :func:`~sunpy.map.header_helper.meta_keywords` can be passed to the :func:`~sunpy.map.header_helper.make_fitswcs_header` and will then populate the returned MetaDict header.
Furthermore, the following observation keywords can be passed to the `~sunpy.map.header_helper.make_fitswcs_header` function: ``observatory``, ``instrument``, ``telescope``, ``wavelength``, ``exposure``.

An example of creating a header with these additional keywords:

.. code-block:: python

    >>> header = make_fitswcs_header(data, coord,
    ...                                        reference_pixel = u.Quantity([5, 5]*u.pixel),
    ...                                        scale = u.Quantity([2, 2] *u.arcsec/u.pixel),
    ...                                        telescope = 'Test case', instrument = 'UV detector',
    ...                                        wavelength = 1000*u.angstrom)
    >>> for key, value in header.items():
    ...     print(f"{key}: {value}")
    wcsaxes: 2
    crpix1: 6.0
    crpix2: 6.0
    cdelt1: 2.0
    cdelt2: 2.0
    cunit1: arcsec
    cunit2: arcsec
    ctype1: HPLN-TAN
    ctype2: HPLT-TAN
    crval1: 0.0
    crval2: 0.0
    lonpole: 180.0
    latpole: 0.0
    mjdref: 0.0
    date-obs: 2013-10-28T00:00:00.000
    rsun_ref: 695700000.0
    dsun_obs: 148644585949.49
    hgln_obs: 0.0
    hglt_obs: 4.7711570596394
    instrume: UV detector
    telescop: Test case
    wavelnth: 1000.0
    waveunit: Angstrom
    naxis: 2
    naxis1: 10
    naxis2: 10
    pc1_1: 1.0
    pc1_2: -0.0
    pc2_1: 0.0
    pc2_2: 1.0
    rsun_obs: 965.3829548285768

From these header MetaDict's that are generated, we can now create a custom map:

.. code-block:: python

    >>> my_map = sunpy.map.Map(data, header)

.. _Meta Keywords Table:

.. list-table:: Meta Keywords
   :widths: 7 30
   :header-rows: 1

   * - Keyword
     - Description
   * - cunit1
     - Units of the coordinate increments along naxis1 e.g. arcsec (required)
   * - cunit2
     - Units of the coordinate increments along naxis2 e.g. arcsec (required)
   * - crval1
     - Coordinate value at reference point on naxis1 (required)
   * - crval2
     - Coordinate value at reference point on naxis2 (required)
   * - cdelt1
     - Spatial scale of pixels for naxis1, i.e. coordinate increment at reference point
   * - cdelt2
     - Spatial scale of pixels for naxis2, i.e. coordinate increment at reference point
   * - crpix1
     - Pixel coordinate at reference point naxis1
   * - crpix2
     - Pixel coordinate at reference point naxis2
   * - ctype1
     - Coordinate type projection along naxis1 of data e.g. HPLT-TAN
   * - ctype2
     - Coordinate type projection along naxis2 of data e.g. HPLN-TAN
   * - hgln_obs
     - Heliographic longitude of observation
   * - hglt_obs
     - Heliographic latitude of observation
   * - dsun_obs
     - distance to Sun from observation in metres
   * - rsun_obs
     - radius of Sun in meters from observation
   * - dateobs
     - date of observation e.g. 2013-10-28 00:00
   * - date_obs
     - date of observation e.g. 2013-10-28 00:00
   * - rsun_ref
     - reference radius of Sun in meters
   * - solar_r
     - radius of Sun in meters from observation
   * - radius
     - radius of Sun in meters from observation
   * - crln_obs
     - Carrington longitude of observation
   * - crlt_obs
     - Heliographic latitude of observation
   * - solar_b0
     - Solar B0 angle
   * - detector
     - name of detector e.g. AIA
   * - exptime
     - exposure time of observation, in seconds e.g 2
   * - instrume
     - name of instrument
   * - wavelnth
     - wavelength of observation
   * - waveunit
     - unit for which observation is taken e.g. angstom
   * - obsrvtry
     - name of observatory of observation
   * - telescop
     - name of telescope of observation
   * - lvl_num
     - FITS processing level
   * - crota2
     - Rotation of the horizontal and vertical axes in degrees
   * - PC1_1
     - Matrix element PCi_j describing the rotation required to align solar North with the top of the image.
   * - PC1_2
     - Matrix element PCi_j describing the rotation required to align solar North with the top of the image.
   * - PC2_1
     - Matrix element PCi_j describing the rotation required to align solar North with the top of the image.
   * - PC2_2
     - Matrix element PCi_j describing the rotation required to align solar North with the top of the image.
   * - CD1_1
     - Matrix element CDi_j describing the rotation required to align solar North with the top of the image.
   * - CD1_2
     - Matrix element CDi_j describing the rotation required to align solar North with the top of the image.
   * - CD2_1
     - Matrix element CDi_j describing the rotation required to align solar North with the top of the image.
   * - CD2_2
     - Matrix element CDi_j describing the rotation required to align solar North with the top of the image.
