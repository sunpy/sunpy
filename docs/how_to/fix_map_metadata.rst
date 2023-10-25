.. _sunpy-how-to-fix-map-metadata:

********************
Fix Metadata for map
********************

There will be times where you will come across either a FITS files with incorrect metadata, or metadata that is not understood by ``sunpy``.
To load these files in to `sunpy.map.Map` you will need to fix the metadata before hand.

This will heavily depend on what metadata is incorrect, for example if the units in the FITS header are incorrect, (which will raise a metadata warning) you can fix this by changing the ``cunit1`` and ``cunit2`` keywords to the correct units:

.. code-block:: python

    >>> from astropy.io import fits

    >>> from sunpy.map import Map
    >>> import sunpy.data.sample

    >>> filepath = sunpy.data.sample.AIA_171_IMAGE  # doctest: +REMOTE_DATA
    >>> data, header = fits.getdata(filepath, header=True)  # doctest: +REMOTE_DATA
    >>> # Note that it is case insensitive for the keys
    >>> header['cunit1'] = 'arcsec'  # doctest: +REMOTE_DATA
    >>> header['cunit2'] = 'arcsec'  # doctest: +REMOTE_DATA
    >>> updated_map = Map(data, header)  # doctest: +REMOTE_DATA

This would apply for the observer_coordinate, wavelength, exposure time, etc.
If you are unsure what the correct metadata should be, you can check the `FITS standard <https://fits.gsfc.nasa.gov/fits_standard.html>`__.
