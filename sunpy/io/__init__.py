"""File input and output functions"""
from __future__ import absolute_import
from sunpy.io import fits, jp2

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
    """Reads the header from a given file"""
    for extension, reader in _known_formats.items():
        if filepath.endswith(extension):
            return reader.get_header(filepath)
        
    reader = detect_filetype(filepath)
    return reader.get_header(filepath)  

def detect_filetype(filepath):
    """Attempts to determine the type of data contained in a file"""
    import re
    
    # Open file and read in first two lines
    with open(filepath) as fp:
        line1 = fp.readline()
        line2 = fp.readline() 
    
    # FITS
    #
    # Checks for "KEY_WORD  =" at beginning of file
    match = re.match(r"[A-Z0-9_]{0,8} *=", line1)
    
    if match is not None and len(match.string) == 9:
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
