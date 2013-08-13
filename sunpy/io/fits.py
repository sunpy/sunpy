"""
FITS File Reader

Notes
-----
FITS
    [1] FITS files allow comments to be attached to every value in the header.
    This is implemented in this module as a KEYCOMMENTS dictionary in the 
    sunpy header. To add a comment to the file on write, add a comment to this
    dictionary with the same name as a key in the header (upcased).

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
import copy
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits

from sunpy.io.header import FileHeader

__all__ = ['read', 'get_header', 'write']

__author__ = "Keith Hughitt, Stuart Mumford"
__email__ = "keith.hughitt@nasa.gov"

def read(filepath):
    """
    Read a fits file
    
    Parameters
    ----------
    filepath : string
        The fits file to be read
        
    Returns
    -------
    pairs : list
        A list of (data, header) tuples
    
    Notes
    -----
    This routine reads all the HDU's in a fits file and returns a list of the 
    data and a FileHeader instance for each one.
    Also all comments in the original file are concatenated into a single
    'comment' key in the returned FileHeader.
    """
    hdulist = pyfits.open(filepath)
    try:
        hdulist.verify('silentfix')
        
        headers = get_header(hdulist)
        pairs = []
        for hdu,header in zip(hdulist, headers):
            pairs.append((hdu.data, header))
    finally:
        hdulist.close()

    return pairs

def get_header(afile):
    """
    Read a fits file and return just the headers for all HDU's
    
    Parameters
    ----------
    afile : string or pyfits.HDUList
        The file to be read, or HDUList to process
    
    Returns
    -------
    headers : list
        A list of FileHeader headers
    """
    if isinstance(afile,pyfits.HDUList):
        hdulist = afile
        close = False
    else:
        hdulist = pyfits.open(afile)
        hdulist.verify('silentfix')
        close=True
        
    try:
        headers= []
        for hdu in hdulist:
            try:
                comment = "".join(hdu.header['COMMENT']).strip()
            except KeyError:
                comment = ""
            try:
                history = "".join(hdu.header['HISTORY']).strip()
            except KeyError:
                history = ""
            
            header = FileHeader(hdu.header)
            header['COMMENT'] = comment
            header['HISTORY'] = history
            
            #Strip out KEYCOMMENTS to a dict, the hard way
            keydict = {}
            for card in hdu.header.cards:
                if card.comment != '':
                 keydict.update({card.keyword:card.comment})
            header['KEYCOMMENTS'] = keydict
            
            headers.append(header)
    finally:
        if close:
            hdulist.close()
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
    #Copy header so the one in memory is left alone while changing it for write
    header = header.copy()
    #The comments need to be added to the header seperately from the normal
    # kwargs. Find and deal with them:
    fits_header = pyfits.Header()
    # Check Header
    key_comments = header.pop('KEYCOMMENTS', False)

    for k,v in header.items():
        if isinstance(v, pyfits.header._HeaderCommentaryCards):
            if k is 'comments':
                fits_header.add_comments(str(v))
            elif k in 'history':
                fits_header.add_history(str(v))
            else:
                fits_header.append(pyfits.Card(k, str(v)))
        else:
            fits_header.append(pyfits.Card(k,v))

    
    if isinstance(key_comments, dict):
            for k,v in key_comments.items():
                fits_header.comments[k] = v
    elif key_comments:
        raise TypeError("KEYCOMMENTS must be a dictionary")
    
    fitskwargs = {'output_verify':'fix'}
    fitskwargs.update(kwargs)
    pyfits.writeto(os.path.expanduser(fname), data, header=fits_header, **fitskwargs)