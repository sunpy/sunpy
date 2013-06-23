from __future__ import absolute_import

from sunpy.io import fits, jp2

__all__ = ['read_file', 'read_file_header', 'write_file', 'detect_filetype']

# File formats supported by SunPy
_known_formats = {
    ('fts', 'fits'): fits,
    ('jp2', 'j2k', 'jpc', 'jpt'): jp2
}

def read_file(filepath):
    """Determines the filetype and reads in the file"""
    for extension, reader in _known_formats.items():
        if filepath.endswith(extension):
            return reader.read(filepath)

    # If filetype is not apparent from extension, attempt to detect
    reader = detect_filetype(filepath)    
    return reader.read(filepath)

def read_file_header(filepath):
    """
    Reads the header from a given file
    
    This should always return a instance of io.header.FileHeader
    """
    for extension, reader in _known_formats.items():
        if filepath.endswith(extension):
            return reader.get_header(filepath)
        
    reader = detect_filetype(filepath)
    return reader.get_header(filepath)  

def write_file(fname, data, header, filetype='auto', **kwargs):
    """
    Write a file from a data & header pair using one of the defined file types
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
    
def detect_filetype(filepath):
    """Attempts to determine the type of data contained in a file"""
    import re
    
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
    raise UnrecognizedFileTypeError("The requested filetype is not currently "
                                    "supported by SunPy.")

class UnrecognizedFileTypeError(IOError):
    """Exception to raise when an unknown file type is encountered"""
    pass

class InvalidJPEG2000FileExtension(IOError):
    """Exception to raise when JPEG 2000 is detected but contains an invalid
    file extension (and is thus unreadable by j2k_to_image"""
    pass
