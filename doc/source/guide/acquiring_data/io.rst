------------------------
Opening Files with SunPy
------------------------

SunPy has wrapped several libraries in order to make input and output as painless as possible. 
Below is a brief tour of the IO libraries wrapped for SunPy. 

1. Wrapped Libraries
====================

Currently three IO libraries are wrapped for easy use with SunPy.

1. Astropy's Fits 
2. Glymur's JPEG 2000 C library wrapper
3. ANA C library

Part of SunPys convenience is that these file types have general purpose high-level IO functions.::

	import sunpy.io as io
	io.read_file(filename)
	io.read_file_header(filename)
	io.write_file(filename)

These are designed to take a filename and breakdown the name in order to automatically calculate the call required to either read or write this file. 
SunPy has a list of known file types which is compared against when using the high-level IO functions.
When reading a data file, the function will return a list of (data, header) pairs depending on how HDUs exist in the file. 
It important to remember this.

Further, you can force the filetype from this interface, like so::

	io.read_file(filename,filetype= 'filetype of your choice')

Valid values for filetype are 'fits', 'jp2' and 'ana'
This will work for the three high-level IO functions.

Full documentation for compatible files is located here :ref:`iofits`, :ref:`iojp2` and :ref:`ioana`.
