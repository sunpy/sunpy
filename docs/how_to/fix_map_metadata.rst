.. _sunpy-how-to-fix-map-metadata:

*************************
Fixing incorrect metadata
*************************

There will be times where you will come across a FITS files with either incorrect, missing or unparsable metadata and reading these files into `~sunpy.map.Map` will cause an error.
Therefore, to load these files into a `~sunpy.map.Map`, you will need to correct the metadata beforehand.

In the example below, the units in the FITS header, as controlled by the ``CUNIT1`` and ``CUNIT2`` keywords, are incorrect.
Before loading the file into a `~sunpy.map.Map`, we will correct these keywords to have the correct units.

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

This applies for any keyword in the `FITS standard <https://fits.gsfc.nasa.gov/fits_standard.html>`__.
