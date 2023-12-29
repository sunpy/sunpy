.. _sunpy-how-to-fix-map-metadata:

*************************
Fixing incorrect metadata
*************************

There will be times where you will come across a FITS files with either incorrect, missing or unparsable metadata and reading these files into `~sunpy.map.Map` will cause an error.
Therefore, to load these files into a `~sunpy.map.Map`, you will need to correct the metadata before hand.

This will heavily depend on what metadata requires changing, unforutanlly it is near impossible to have a guide for every possible keyword.
Below we have an example where the units in the FITS header are incorrect, this is controlled by the ``cunit1`` and ``cunit2`` keywords in a FITS header.
So we will change those values and then load the file into a `~sunpy.map.Map`.

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

This applies for any FITS standard keyword, which you can find in the `FITS standard <https://fits.gsfc.nasa.gov/fits_standard.html>`__.
