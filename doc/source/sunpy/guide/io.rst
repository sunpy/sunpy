============================
Your Guide to Opening Files!
============================

Some of have data. Now, this data is most likely on a PC somewhere in the world (unless you have non-digital data which in case you should stop reading.)
SunPy has several inbuilt or wrapped functions calls around libraries that open a variety of files. 
Let me take you on a journey through the binary forests. 


Wrapped Libraries
=================

Currently three IO libraries are wrapped for easy use with SunPy.

1. AstroPy's Fits 
2. Glymur's JPEG 2000 C library wrapper
3. ANA C library

Luckily for everyone, SunPy has general purpose high-level IO functions.::

	import sunpy.io as io
	io.file_tiils.read_file(filename)
	io.file_tools.read_file_header(filename)
	io.file_tools.write_file(filename)

These are designed to take a filename and breakdown the name in order to automatically calculate the calls needed to either read or write this file. 
When reading or writing a data file, the function will return a (data, header) pair. 
It important to remember this.
If for some reason these calls fail, which is possible to do unknown file extensions or  unsupported file types, create an issue and someone will happily get back to you.

Sometimes, you require fine control over your data files. 
Come on and follow me.

Fits files
----------

SunPy's Fits reading ability comes directly from AstroPy as they swallowed whole PyFits a while back now. 
Long Live AstroPy.

Currently four functions exist under the io.fits namespace.::
	
	import sunpy.io.fits as fits
	fits.read(filename)
	fits.write(filename)
	fits.get_header(filename)
	fits.extract_waveunit(filename)



Full documenation is located `here <https://sunpy.readthedocs.org/en/latest/reference/generated/sunpy.io.fits.html>`.

JPEG files
-----------

SunPy's JPEG reading ability comes from wrappying Glymur which itself wraps a OpenJPEG 2000 C library. Lots of wrapping here. Christmas come early.

Currently four functions exist under the io.fits namespace.::

	import sunpy.jp2 as jp2
	io.file_tiils.read_file(filename)
	io.file_tools.read_file_header(filename)
	io.file_tools.write_file(filename)

Full documenation is located `here <https://sunpy.readthedocs.org/en/latest/reference/generated/sunpy.io.jp2.html>`.

fz files
----------

SunPy's ANA (fz) reading ability comes from wrapping a C library designed to read and write ANA compressed files. 
It is common to find fits files compressed this way. The fits wrapper already takes care of these files. 

Currently four functions exist under the io.fits namespace.::
	import sunpy.io.ana as ana
	io.file_tiils.read_file(filename)
	io.file_tools.read_file_header(filename)
	io.file_tools.write_file(filename)

Full documenation is located `here <https://sunpy.readthedocs.org/en/latest/reference/generated/sunpy.io.ana.html>`.