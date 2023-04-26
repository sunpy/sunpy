.. _how-to-create-a-map:

How to create a `Map <sunpy.map.GenericMap>`
============================================

The SunPy Map factory accepts a wide variety of inputs for creating maps

* Preloaded tuples of (data, header) pairs

>>> mymap = sunpy.map.Map((data, header))   # doctest: +SKIP

headers are some base of `dict` or `collections.OrderedDict`, including
`sunpy.io.header.FileHeader` or `sunpy.util.metadata.MetaDict` classes.

* data, header pairs, not in tuples

>>> mymap = sunpy.map.Map(data, header)   # doctest: +SKIP

* data, wcs object, in tuple

>>> from astropy.wcs import WCS
>>> wcs = WCS(sunpy.data.sample.AIA_171_ROLL_IMAGE)  # doctest: +SKIP
>>> data = fits.getdata(sunpy.data.sample.AIA_171_ROLL_IMAGE)  # doctest: +SKIP
>>> mymap = sunpy.map.Map((data, wcs))  # doctest: +SKIP

* data, wcs object, not in tuple

>>> from astropy.wcs import WCS
>>> wcs = WCS(sunpy.data.sample.AIA_171_ROLL_IMAGE)  # doctest: +SKIP
>>> data = fits.getdata(sunpy.data.sample.AIA_171_ROLL_IMAGE)  # doctest: +SKIP
>>> mymap = sunpy.map.Map(data, wcs)  # doctest: +SKIP

* File names

>>> mymap = sunpy.map.Map('file1.fits')   # doctest: +SKIP

* All fits files in a directory by giving a directory

>>> mymap = sunpy.map.Map('local_dir/sub_dir')   # doctest: +SKIP

* A filesystem path expressed as a `pathlib.Path`

>>> import pathlib
>>> mymap = sunpy.map.Map(pathlib.Path('file1.fits'))  # doctest: +SKIP
>>> sub_dir = pathlib.Path('local_dir/sub_dir')
>>> mymap = sunpy.map.Map(sub_dir)   # doctest: +SKIP
>>> mymap = sunpy.map.Map(sub_dir / 'file3.fits')   # doctest: +SKIP

* Some regex globs

>>> mymap = sunpy.map.Map('eit_*.fits')   # doctest: +SKIP

* URLs

>>> mymap = sunpy.map.Map(url_str)   # doctest: +SKIP

* DatabaseEntry

>>> mymap = sunpy.map.Map(db_result)   # doctest: +SKIP

* Lists of any of the above

>>> mymap = sunpy.map.Map(['file1.fits', 'file2.fits', 'file3.fits', 'directory1/'])  # doctest: +SKIP

* Any mixture of the above not in a list

>>> mymap = sunpy.map.Map(((data, header), data2, header2, 'file1.fits', url_str, 'eit_*.fits'))  # doctest: +SKIP
