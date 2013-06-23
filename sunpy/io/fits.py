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

import os

import pyfits

from sunpy.io.header import FileHeader

__all__ = ['read', 'get_header']

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

def read(filepath):
    """Reads in the file at the specified location"""
    hdulist = pyfits.open(filepath)
    hdulist.verify('silentfix')
    
    pairs = []
    for hdu in hdulist:
        fits_comment = hdulist[0].header.get_comment()
        
        # PyFITS 2.x
        if len(fits_comment) > 0 and isinstance(fits_comment[0], basestring):
            comments = [val for val in fits_comment]       
        else:
            # PyFITS 3.x
            comments = [card.value for card in fits_comment]
            
        comment = "".join(comments).strip()
        header = FileHeader(hdu.header)
        header['comment'] = comment
        pairs.append((hdu.data, header))

    return pairs

def get_header(filepath):
    """Returns a list of headers for all HDUs"""
    hdulist = pyfits.open(filepath)
    hdulist.verify('silentfix')
    headers= []
    for hdu in hdulist:
        comment = "".join(hdulist[0].header.get_comment()).strip()
        header = FileHeader(hdulist[0].header)
        header['comment'] = comment
        headers.append(header)
    return headers

def write(fname, data, header, **kwargs):
    """
    Take a data header pair and write a fits file
    
    Parameters
    ----------
    fname: str
        File name, with extension
        
    data: ndarray
        n-dimensional data array
    
    header: dict
        A header dictionary
    """
    #The comments need to be added to the header seperately from the normal
    # kwargs. Find and deal with them:
    cards = []
    comments = []
    # Check Header
    for k,v in header.items():
        if isinstance(v, pyfits.header._HeaderCommentaryCards):
            comments.append(v)
        else:
            cards.append(pyfits.core.Card(k,v))
    cards = pyfits.core.Header(cards)
    if len(comments) > 0:
        cards.add_comment(comments)
    
    fitskwargs = {'output_verify':'fix'}
    fitskwargs.update(kwargs)
    pyfits.writeto(os.path.expanduser(fname), data, header=cards, **fitskwargs)