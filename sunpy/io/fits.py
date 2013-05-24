"""
FITS File Reader

Notes
-----
PyFITS
    [1] Due to the way PyFITS works with images the header dictionary may
    differ depending on whether is accessed before or after the fits[0].data
    is requested. If the header is read before the data then the original
    header will be returned. If the header is read after the data has been
    accessed then the data will have been scaled and a modified header
    reflecting these changes will be returned: BITPIX may differ and
    BSCALE and B_ZERO may be dropped in the modified version.
    
    [2] The verify('fix') call attempts to handle violations of the FITS
    standard. For example, nan values will be converted to "nan" strings.
    Attempting to cast a pyfits header to a dictionary while it contains
    invalid header tags will result in an error so verifying it early on
    makes the header easier to work with later.

References
----------
| http://stackoverflow.com/questions/456672/class-factory-in-python
| http://stsdas.stsci.edu/download/wikidocs/The_PyFITS_Handbook.pdf

"""
from __future__ import absolute_import

import pyfits

from sunpy.io.header import FileHeader

__all__ = ['read', 'get_header']

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

class FITSHeader(FileHeader):
    """
    FITSHeader(header)

    A dictionary-like class for working with FITS  headers. Upper cases all
    keywords to conform to FITSstandard

    Parameters
    ----------
    header : pyfits.core.Header, dict
        Header tags associated with the data

    """
    def __init__(self, adict):
        """Creates a new MapHeader instance"""
        # Store all keys as upper-case to allow for case-insensitive indexing
        adict = dict((k.upper(), v) for k, v in adict.items())
        #TODO: This might not work!
        FileHeader.__init__(self, adict)
        
    def __contains__(self, key):
        """Overide __contains__"""
        return dict.__contains__(self, key.upper())

    def __getitem__(self, key):
        """Overide [] indexing"""
        return dict.__getitem__(self, key.upper())

    def __setitem__(self, key, value):
        """Overide [] indexing"""
        return dict.__setitem__(self, key.upper(), value)
    
    def as_pyfits_header(self):
        """Returns a PyFITS header instance of the header"""
        cards = [pyfits.core.Card(k, v) for k, v in self.items()]
        return pyfits.core.Header(cards)

    def copy(self):
        """Overide copy operator"""
        return type(self)(dict.copy(self))

    def get(self, key, default=None):
        """Overide .get() indexing"""
        return dict.get(self, key.upper(), default)

    def has_key(self, key):
        """Overide .has_key() to perform case-insensitively"""
        return key.upper() in self

    def pop(self, key, default=None):
        """Overide .pop() to perform case-insensitively"""
        return dict.pop(self, key.upper(), default)

    def update(self, d2):
        """Overide .update() to perform case-insensitively"""
        return dict.update(self, dict((k.upper(), v) for k, v in d2.items()))

    def setdefault(self, key, default=None):
        """Overide .setdefault() to perform case-insensitively"""
        return dict.setdefault(self, key.upper(), default)

def read(filepath):
    """Reads in the file at the specified location"""
    hdulist = pyfits.open(filepath)
    hdulist.verify('silentfix')
    
    fits_comment = hdulist[0].header.get_comment()
    
    # PyFITS 2.x
    if len(fits_comment) > 0 and isinstance(fits_comment[0], basestring):
        comments = [val for val in fits_comment]       
    else:
        # PyFITS 3.x
        comments = [card.value for card in fits_comment]
        
    comment = "".join(comments).strip()
    header = FITSHeader(hdulist[0].header)
    header['comment'] = comment

    return hdulist[0].data, header

def get_header(filepath):
    """Returns the header for a given file"""
    hdulist = pyfits.open(filepath)
    hdulist.verify('silentfix')
    
    comment = "".join(hdulist[0].header.get_comment()).strip()
    header = FITSHeader(hdulist[0].header)
    header['comment'] = comment
            
    return header
