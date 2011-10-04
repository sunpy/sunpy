from __future__ import absolute_import
"""
File input and output functions
"""
from sunpy.map.io import fits

# File formats supported by SunPy
_known_formats = {
    ('fts', 'fits'): fits
}

def read_file(filepath):
    """Determines the filetype and reads in the file"""
    for extension, reader in _known_formats.items():
        if filepath.endswith(extension):
            return reader.read(filepath)
