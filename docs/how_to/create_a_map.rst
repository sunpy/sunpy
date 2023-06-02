.. _sunpy-how-to-create-a-map:

*************************
How to create a sunpy Map
*************************

One of the primary goals of the Map interface is to make it as easy as possible to create a Map.
As such, you can pass many different kinds of inputs to Map.
These are listed below.

File name
=========

If you have a FITS file, this is the easiest and recommended way to create a Map.
This can be either a string or a `~pathlib.Path`.

.. code-block:: python

    >>> import pathlib

    >>> import sunpy.map
    >>> import sunpy.data.sample

    >>> my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
    >>> my_map = sunpy.map.Map('file.fits')  # doctest: +SKIP
    >>> my_map = sunpy.map.Map(pathlib.Path('file.fits'))  # doctest: +SKIP
    >>> sub_dir = pathlib.Path('local_dir/sub_dir')
    >>> my_map = sunpy.map.Map(sub_dir / 'another_file.fits')   # doctest: +SKIP

Directory containing FITS files
===============================

If there is more than one FITS file in the directory, this will return a list of Map objects.

.. code-block:: python

    >>> my_maps = sunpy.map.Map('local_dir/sub_dir')   # doctest: +SKIP
    >>> my_maps = sunpy.map.Map(sub_dir)   # doctest: +SKIP

Array and `astropy.io.fits.Header`
==================================

If needed, this way can be used to modify the header before passing it to `~sunpy.map.Map`.

.. code-block:: python

    >>> import astropy.io.fits

    >>> with astropy.io.fits.open(sunpy.data.sample.AIA_171_IMAGE) as hdul:
    ...     data = hdul[1].data
    ...     header = hdul[1].header  # doctest: +REMOTE_DATA
    >>> my_map = sunpy.map.Map(data, header)  # doctest: +REMOTE_DATA

These data header pairs can also be passed as a `tuple`,

.. code-block:: python

    >>> my_map = sunpy.map.Map((data, header))  # doctest: +REMOTE_DATA

Data array and a `~sunpy.util.metadata.MetaDict` object
=======================================================

This includes any base class of `~sunpy.util.metadata.MetaDict`, including `dict` or `collections.OrderedDict`.

.. code-block:: python

    >>> import sunpy.util.metadata

    >>> meta = sunpy.util.metadata.MetaDict(header)  # doctest: +REMOTE_DATA
    >>> my_map = sunpy.map.Map(data, meta)   # doctest: +REMOTE_DATA

Data array and an `astropy.wcs.WCS` object
==========================================

.. code-block:: python

    >>> import astropy.wcs

    >>> wcs = astropy.wcs.WCS(header=header)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
    >>> my_map = sunpy.map.Map(data, wcs)  # doctest: +REMOTE_DATA

Glob patterns
=============

If the glob pattern matches more than one FITS file, this will return a list of Map objects.

.. code-block:: python

    >>> my_map = sunpy.map.Map('eit_*.fits')   # doctest: +SKIP

URL
===

.. code-block:: python

    >>> sample_data_url = 'http://data.sunpy.org/sunpy/v1/AIA20110607_063302_0171_lowres.fits'
    >>> my_map = sunpy.map.Map(sample_data_url)  # doctest: +REMOTE_DATA

Combinations of any of the above
================================

These can either be in a list or as separate arguments.
As with the case of a directory or glob pattern, this will return multiple Map objects.

.. code-block:: python

    >>> my_map = sunpy.map.Map(['file1.fits', 'file2.fits', 'file3.fits', 'directory1/'])  # doctest: +SKIP
    >>> my_map = sunpy.map.Map((data, header), data, meta, 'file1.fits', sample_data_url, 'eit_*.fits')  # doctest: +SKIP
