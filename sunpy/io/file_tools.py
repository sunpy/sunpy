from __future__ import absolute_import

import re

# File formats supported by SunPy
try:
    from sunpy.io import fits
    _fits_imported = True
except ImportError:
    _fits_imported = False

try:
    from sunpy.io import jp2
    _jp2_imported = True
except ImportError:
    _jp2_imported = False

__all__ = ['read_file', 'read_file_header', 'write_file']

_known_formats = {
    ('fts', 'fits'): 'fits',
    ('jp2', 'j2k', 'jpc', 'jpt'): 'jp2'
}
def read_file(filepath, **kwargs):
    """
    Automatically determine the filetype and read the file
    
    Parameters
    ----------
    filepath : string
        The file to be read
    
    Returns
    -------
    pairs : list
        A list of (data, header) tuples.
    """
    # If filetype is not apparent from extension, attempt to detect
    ftype = _detect_filetype(filepath)
    reader = _get_reader(ftype)
    return reader.read(filepath, **kwargs)

def read_file_header(filepath, **kwargs):
    """
    Reads the header from a given file
    
    This should always return a instance of io.header.FileHeader
    
    Parameters
    ----------
    
    filepath :  string
        The file from which the header is to be read.
    
    Returns
    -------
    
    headers : list
        A list of headers
    """
    ftype = _detect_filetype(filepath)
    reader = _get_reader(ftype)
    
    return reader.get_header(filepath, **kwargs)  

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
    
    Other keyword arguments will be passes to the writer function used.
    
    This routine currently only supports saving a single HDU.
    """
    if filetype == 'auto':
        ftype = _detect_filetype(fname, write=True)
        reader = _get_reader(ftype)
    
    else:
        for extension, ftype in _known_formats.items():
            if filetype in extension:
                reader = _get_reader(ftype)
                return reader.write(fname, data, header, **kwargs)
            
    #Nothing has matched, panic
    raise ValueError("This filetype is not supported")
    
def _get_reader(ftype):
    if ftype is 'fits':
        if _fits_imported:
            return fits
        else:
            raise ImportError(
            "FITS reader can not be imported, do you have astropy.io.fits?")
        
    if ftype is 'jp2':
        if _jp2_imported:
            return jp2
        else:
            #This one will not be used until we have a glymur dep
            raise ImportError(
            "jp2 reader can not be imported, do you have openjpeg installed?")
    
def _detect_filetype(filepath, write=False):
    ftype = None
    for extension, reader in _known_formats.items():
        if filepath.endswith(extension):
            ftype = reader
    
    if not ftype and not write:
        ftype = _detect_filetype_binary(filepath)

    return ftype
    
def _detect_filetype_binary(filepath):
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
            # j2k_to_image requires a valid extension
            raise InvalidJPEG2000FileExtension
            #return 'jp2'

class UnrecognizedFileTypeError(IOError):
    """Exception to raise when an unknown file type is encountered"""
    pass

class InvalidJPEG2000FileExtension(IOError):
    """Exception to raise when JPEG 2000 is detected but contains an invalid
    file extension (and is thus unreadable by j2k_to_image"""
    pass
