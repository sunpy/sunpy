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

from sunpy.map.header import MapHeader

__all__ = ['read', 'get_header']

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

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
    header = MapHeader(hdulist[0].header)
    header['comment'] = comment

    return hdulist[0].data, header

def get_header(filepath):
    """Returns the header for a given file"""
    hdulist = pyfits.open(filepath)
    hdulist.verify('silentfix')
    
    comment = "".join(hdulist[0].header.get_comment()).strip()
    header = MapHeader(hdulist[0].header)
    header['comment'] = comment
            
    return header
