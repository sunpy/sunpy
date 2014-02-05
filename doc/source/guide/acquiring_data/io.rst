------------------------
Opening Files with SunPy
------------------------

SunPy has wrapped several libraries in order to make input and output as painless as possible. 
So let me take you on a journey through the binary forests using the chainsaw of SunPy. 

1. Wrapped Libraries
====================

Currently three IO libraries are wrapped for easy use with SunPy.

1. AstroPy's Fits 
2. Glymur's JPEG 2000 C library wrapper
3. ANA C library

Luckily for everyone, SunPy has general purpose high-level IO functions.::

	import sunpy.io as io
	io.file_tools.read_file(filename)
	io.file_tools.read_file_header(filename)
	io.file_tools.write_file(filename)

These are designed to take a filename and breakdown the name in order to automatically calculate the call required to either read or write this file. 
SunPy has a "database" of file extensions which is compared against when using the high-level IO functions.
When reading or writing a data file, the function will return a list of (data, header) pairs depending on how HDUs exist in the file. 
It important to remember this.
If for some reason these calls fail, which is possible to do unknown file extensions or unsupported file types, create an issue and someone will get back to you.

Further, you can force the filetype from this interface, like so::

	io.fileTools.read_file(filename,filetype= 'filetype of your choice')

This will work for the three high-level IO functions.

Sometimes, you require fine control over your data files. 
Come on and follow me deeper.

1.1 Fits files
--------------

SunPy's Fits reading ability comes directly from AstroPy, as they swallowed PyFits whole a while back now. 

Currently four functions exist under the io.fits namespace.::
	
	import sunpy.io.fits as fits
	fits.read(filename)
	fits.write(filename, data, header)
	fits.get_header(filename)
	fits.extract_waveunit(header)

So, the functions here are basic. You can read/write to/from fit files but also get only the header as well.
The final function is used to figure out the wave unit from a fits header and so you should pass a header object into this function.

Full documentation is located `here <https://sunpy.readthedocs.org/en/latest/reference/generated/sunpy.io.fits.html>`.

1.2 JPEG 2000 files
-------------------

SunPy's JPEG reading ability comes from wrappying Glymur which itself wraps a OpenJPEG 2000 C library. Lots of wrapping here. Christmas come early.

Currently four functions exist under the io.jp2 namespace.::

	import sunpy.jp2 as jp2
	jp2.read_file(filename)
	jp2.read_file_header(filename)
	jp2.write_file(filename)

Full documentation is located `here <https://sunpy.readthedocs.org/en/latest/reference/generated/sunpy.io.jp2.html>`.

1.3 fz files
------------

SunPy's ANA (fz) reading ability comes from wrapping a C library designed to read and write ANA compressed files. 
It is common to find fits files compressed this way. The fits wrapper already takes care of these files. 

Currently four functions exist under the io.ana namespace.::

	import sunpy.io.ana as ana
	ana.read_file(filename)
	ana.read_file_header(filename)
	ana.write_file(filename)

Full documentation is located `here <https://sunpy.readthedocs.org/en/latest/reference/generated/sunpy.io.ana.html>`.

1.4 Other files
---------------

For these, you are on your own. 
However, Python has many IO libraries and there will be the ability to read them in.
You can join the IRC channel for support on freenode (#sunpy) or the mailing list (link) and the github issues section (link).