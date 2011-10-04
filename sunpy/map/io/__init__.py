from __future__ import absolute_import
"""
File input and output functions
"""
from sunpy.map.io import fits, jp2

# File formats supported by SunPy
_known_formats = {
    ('fts', 'fits'): fits,
    ('jp2'): jp2
}

def read_file(filepath):
    """Determines the filetype and reads in the file"""
    for extension, reader in _known_formats.items():
        if filepath.endswith(extension):
            return reader.read(filepath)

    # Raise an error if an unsupported filetype is encountered
    raise UnrecognizedFileTypeError

class UnrecognizedFileTypeError(IOError):
    """Exception to raise when an unknown file type is encountered"""
    pass
