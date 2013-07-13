from __future__ import absolute_import

import re
import warnings

# File formats supported by SunPy
# If dependancy is not importable, do not add it to the readers.
_known_formats = {}
try:
    from sunpy.io import fits
    _known_formats.update({('fts', 'fits'): fits})
    _fits_failed = False
except ImportError:
    _fits_failed = True
    warnings.warn(Warning("FITS reader failed to import"))

try:
    from sunpy.io import jp2
    _known_formats.update({('jp2', 'j2k', 'jpc', 'jpt'): jp2})
    _jp2_failed = False
except ImportError:
    _jp2_failed = True
    warnings.warn(Warning("jp2 reader failed to import"))

__all__ = ['read_file', 'read_file_header', 'write_file']

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
    for extension, reader in _known_formats.items():
        if filepath.endswith(extension):
            return reader.read(filepath, **kwargs)

    # If filetype is not apparent from extension, attempt to detect
    reader = _detect_filetype(filepath)    
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
    for extension, reader in _known_formats.items():
        if filepath.endswith(extension):
            return reader.get_header(filepath, **kwargs)
        
    reader = _detect_filetype(filepath)
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
        for extension, reader in _known_formats.items():
            if fname.endswith(extension):
                return reader.write(fname, data, header, **kwargs)
    
    else:
        for extension, reader in _known_formats.items():
            if filetype in extension:
                return reader.write(fname, data, header, **kwargs)
            
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
        return fits
    
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

    # Raise an error if an unsupported filetype is encountered
    if _fits_failed and _jp2_failed:
        raise ImportError(
        """SunPy's FITS reader and JPEG2000 reader failed to import
        FITS reading requres astropy.io and JPEG2000 requires glymur and 
        openjpeg>=1.5""")
    elif _fits_failed:
        raise ImportError(
        "SunPy's FITS reader failed to import FITS reading requres astropy.io")
    elif _jp2_failed:
        raise ImportError(
        """SunPy's JPEG2000 reader failed to import JPEG2000 requires glymur 
        and openjpeg>=1.5""")
    else:
        raise UnrecognizedFileTypeError(
                  "SunPy can not recoginse or does not support this file type")

class UnrecognizedFileTypeError(IOError):
    """Exception to raise when an unknown file type is encountered"""
    pass

class InvalidJPEG2000FileExtension(IOError):
    """Exception to raise when JPEG 2000 is detected but contains an invalid
    file extension (and is thus unreadable by j2k_to_image"""
    pass
