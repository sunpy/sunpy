from __future__ import absolute_import

import re

try:
    from . import fits
except ImportError:
    fits = None

try:
    from . import jp2
except ImportError:
    jp2 = None

try:
    from . import ana
except ImportError:
    ana = None

__all__ = ['read_file', 'read_file_header', 'write_file']

# File formats supported by SunPy
_known_extensions = {
    ('fts', 'fits'): 'fits',
    ('jp2', 'j2k', 'jpc', 'jpt'): 'jp2',
    ('fz', 'f0'): 'ana'
}

#Define a dict which raises a custom error message if the value is None
class Readers(dict):
    def __init__(self, *args):
        dict.__init__(self, *args)

    def __getitem__(self, key):
        val = dict.__getitem__(self, key)
        if val is None:
            raise ReaderError("The Reader sunpy.io.%s is not avalible, please check that you have the required dependancies installed."%key)
        return val

#Map the readers
_readers = Readers({
            'fits':fits,
            'jp2':jp2,
            'ana':ana
})

def read_file(filepath, filetype=None, **kwargs):
    """
    Automatically determine the filetype and read the file

    Parameters
    ----------
    filepath : string
        The file to be read

    filetype: string
        Supported reader or extension to manually specify the filetype.
        Supported readers are ('jp2', 'fits', 'ana')

    Returns
    -------
    pairs : list
        A list of (data, header) tuples.
    """
    if filetype:
        return _readers[filetype].read(filepath, **kwargs)

    for extension, readername in _known_extensions.items():
        if filepath.endswith(extension) or filetype in extension:
            return _readers[readername].read(filepath, **kwargs)

    # If filetype is not apparent from extension, attempt to detect
    readername = _detect_filetype(filepath)
    return _readers[readername].read(filepath, **kwargs)

def read_file_header(filepath, filetype=None, **kwargs):
    """
    Reads the header from a given file

    This should always return a instance of io.header.FileHeader

    Parameters
    ----------

    filepath :  string
        The file from which the header is to be read.

    filetype: string
        Supported reader or extension to manually specify the filetype.
        Supported readers are ('jp2', 'fits')

    Returns
    -------

    headers : list
        A list of headers
    """
    if filetype:
        return _readers[filetype].get_header(filepath, **kwargs)

    for extension, readername in _known_extensions.items():
        if filepath.endswith(extension) or filetype in extension:
            return _readers[readername].get_header(filepath, **kwargs)

    readername = _detect_filetype(filepath)
    return _readers[readername].get_header(filepath, **kwargs)

def write_file(fname, data, header, filetype='auto', **kwargs):
    """
    Write a file from a data & header pair using one of the defined file types.

    Parameters
    ----------
    fname : string
        Filename of file to save

    data : ndarray
        Data to save to a fits file

    header : OrderedDict
        Meta data to save with the data

    filetype : string
        {'auto', 'fits', 'jp2'} Filetype to savem if auto fname extension will
        be detected, else specifiy a supported file extension.
    
    Notes
    -----
    * Other keyword arguments will be passes to the writer function used.
    * This routine currently only supports saving a single HDU.
    """
    if filetype == 'auto':
        for extension, readername in _known_extensions.items():
            if fname.endswith(extension):
                return _readers[readername].write(fname, data, header, **kwargs)

    else:
        for extension, readername in _known_extensions.items():
            if filetype in extension:
                return _readers[readername].write(fname, data, header, **kwargs)

    #Nothing has matched, panic
    raise ValueError("This filetype is not supported" )

def _detect_filetype(filepath):
    """
    Attempts to determine the type of data contained in a file.

    This is only used for reading because it opens the file to check the data.
    """

    # Open file and read in first two lines
    with open(filepath) as fp:
        line1 = fp.readline()
        line2 = fp.readline()
        #Some FITS files do not have line breaks at the end of header cards.
        fp.seek(0)
        first80 = fp.read(80)

    # FITS
    #
    # Checks for "KEY_WORD  =" at beginning of file
    match = re.match(r"[A-Z0-9_]{0,8} *=", first80)

    if match is not None:
        return 'fits'

    # JPEG 2000
    #
    # Checks for one of two signatures found at beginning of all JP2 files.
    # Adapted from ExifTool
    # [1] http://www.sno.phy.queensu.ca/~phil/exiftool/
    # [2] http://www.jpeg.org/public/fcd15444-2.pdf
    # [3] ftp://ftp.remotesensing.org/jpeg2000/fcd15444-1.pdf
    jp2_signatures = ["\x00\x00\x00\x0cjP  \x0d\x0a\x87\x0a",
                      "\x00\x00\x00\x0cjP\x1a\x1a\x0d\x0a\x87\x0a"]

    for sig in jp2_signatures:
        if line1 + line2 == sig:
            return 'jp2'

    # Raise an error if an unsupported filetype is encountered
    raise UnrecognizedFileTypeError("The requested filetype is not currently "
                                    "supported by SunPy.")

class UnrecognizedFileTypeError(IOError):
    """Exception to raise when an unknown file type is encountered"""
    pass

class ReaderError(ImportError):
    """Exception to raise when an unknown file type is encountered"""
    pass

class InvalidJPEG2000FileExtension(IOError):
    pass