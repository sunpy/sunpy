"""
FITS File Reader

Notes
-----
FITS
    [1] FITS files allow comments to be attached to every value in the header.
    This is implemented in this module as a KEYCOMMENTS dictionary in the
    sunpy header. To add a comment to the file on write, add a comment to this
    dictionary with the same name as a key in the header (upcased).

    [2] Due to the way `~astropy.io.fits` works with images the header dictionary may
    differ depending on whether is accessed before or after the fits[0].data
    is requested. If the header is read before the data then the original
    header will be returned. If the header is read after the data has been
    accessed then the data will have been scaled and a modified header
    reflecting these changes will be returned: BITPIX may differ and
    BSCALE and B_ZERO may be dropped in the modified version.

    [3] The verify('fix') call attempts to handle violations of the FITS
    standard. For example, nan values will be converted to "nan" strings.
    Attempting to cast a pyfits header to a dictionary while it contains
    invalid header tags will result in an error so verifying it early on
    makes the header easier to work with later.

References
----------
| https://stackoverflow.com/questions/456672/class-factory-in-python
"""
from __future__ import absolute_import, division, print_function
import os
import re
import sys
import warnings
import traceback
import itertools
import collections

from astropy.io import fits

from sunpy.io.header import FileHeader
from sunpy.extern.six.moves import zip

__all__ = ['read', 'get_header', 'write', 'extract_waveunit']

__author__ = "Keith Hughitt, Stuart Mumford, Simon Liedtke"
__email__ = "keith.hughitt@nasa.gov"

HDPair = collections.namedtuple('HDPair', ['data', 'header'])


def read(filepath, hdus=None, memmap=None, **kwargs):
    """
    Read a fits file

    Parameters
    ----------
    filepath : `str`
        The fits file to be read
    hdu: `int` or iterable
        The HDU indexes to read from the file

    Returns
    -------
    pairs : `list`
        A list of (data, header) tuples

    Notes
    -----
    This routine reads all the HDU's in a fits file and returns a list of the
    data and a FileHeader instance for each one.
    Also all comments in the original file are concatenated into a single
    'comment' key in the returned FileHeader.
    """
    with fits.open(filepath, ignore_blank=True, memmap=memmap) as hdulist:
        if hdus is not None:
            if isinstance(hdus, int):
                hdulist = hdulist[hdus]
            elif isinstance(hdus, collections.Iterable):
                hdulist = [hdulist[i] for i in hdus]

        hdulist.verify('silentfix+warn')

        headers = get_header(hdulist)
        pairs = []

        for i, (hdu, header) in enumerate(zip(hdulist, headers)):
            try:
                pairs.append(HDPair(hdu.data, header))
            except (KeyError, ValueError) as e:
                message = "Error when reading HDU {}. Skipping.\n".format(i)
                for line in traceback.format_tb(sys.exc_info()[2]):
                    message += line
                    message += '\n'
                message += repr(e)
                warnings.warn(message, Warning, stacklevel=2)

    return pairs


def get_header(afile):
    """
    Read a fits file and return just the headers for all HDU's. In each header,
    the key WAVEUNIT denotes the wavelength unit which is used to describe the
    value of the key WAVELNTH.

    Parameters
    ----------
    afile : `str` or fits.HDUList
        The file to be read, or HDUList to process.

    Returns
    -------
    headers : `list`
        A list of FileHeader headers.
    """
    if isinstance(afile, fits.HDUList):
        hdulist = afile
        close = False
    else:
        hdulist = fits.open(afile, ignore_blank=True)
        hdulist.verify('silentfix')
        close = True

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

            # Strip out KEYCOMMENTS to a dict, the hard way
            keydict = {}
            for card in hdu.header.cards:
                if card.comment != '':
                    keydict.update({card.keyword:card.comment})
            header['KEYCOMMENTS'] = keydict
            header['WAVEUNIT'] = extract_waveunit(header)

            headers.append(header)
    finally:
        if close:
            hdulist.close()
    return headers


def write(fname, data, header, **kwargs):
    """
    Take a data header pair and write a FITS file.

    Parameters
    ----------
    fname : `str`
        File name, with extension

    data : `numpy.ndarray`
        n-dimensional data array

    header : `dict`
        A header dictionary
    """
    # Copy header so the one in memory is left alone while changing it for
    # write.
    header = header.copy()

    # The comments need to be added to the header separately from the normal
    # kwargs. Find and deal with them:
    fits_header = fits.Header()
    # Check Header
    key_comments = header.pop('KEYCOMMENTS', False)

    for k, v in header.items():
        if isinstance(v, fits.header._HeaderCommentaryCards):
            if k == 'comments':
                comments = str(v).split('\n')
                for com in comments:
                    fits_header.add_comments(com)
            elif k == 'history':
                hists = str(v).split('\n')
                for hist in hists:
                    fits_header.add_history(hist)
            elif k != '':
                fits_header.append(fits.Card(k, str(v).split('\n')))

        else:
            fits_header.append(fits.Card(k, v))

    if isinstance(key_comments, dict):
        for k, v in key_comments.items():
            fits_header.comments[k] = v
    elif key_comments:
        raise TypeError("KEYCOMMENTS must be a dictionary")

    if isinstance(fname, str):
        fname = os.path.expanduser(fname)

    fitskwargs = {'output_verify': 'fix'}
    fitskwargs.update(kwargs)
    fits.writeto(fname, data, header=fits_header, **fitskwargs)


def extract_waveunit(header):
    """Attempt to read the wavelength unit from a given FITS header.

    Parameters
    ----------
    header : FileHeader
        One :class:`sunpy.io.header.FileHeader` instance which was created by
        reading a FITS file. :func:`sunpy.io.fits.get_header` returns a list of
        such instances.

    Returns
    -------
    waveunit : `str`
        The wavelength unit that could be found or ``None`` otherwise.

    Examples
    --------
    The goal of this function is to return a string that can be used in
    conjunction with the astropy.units module so that the return value can be
    directly passed to ``astropy.units.Unit``::

        >>> import astropy.units
        >>> header = {'WAVEUNIT': 'Angstrom', 'KEYCOMMENTS': {}}
        >>> waveunit = extract_waveunit(header)
        >>> if waveunit is not None:
        ...     unit = astropy.units.Unit(waveunit)

    """
    # algorithm: try the following procedures in the following order and return
    # as soon as a waveunit could be detected
    # 1. read header('WAVEUNIT'). If None, go to step 2.
    # 1.1 -9 -> 'nm'
    # 1.2 -10 -> 'angstrom'
    # 1.3 0 -> go to step 2
    # 1.4 if neither of the above, return the value itself in lowercase
    # 2. parse waveunit_comment
    # 2.1 'in meters' -> 'm'
    # 3. parse wavelnth_comment
    # 3.1 "[$UNIT] ..." -> $UNIT
    # 3.2 "Observed wavelength ($UNIT)" -> $UNIT
    def parse_waveunit_comment(waveunit_comment):
        if waveunit_comment == 'in meters':
            return 'm'

    waveunit_comment = header['KEYCOMMENTS'].get('WAVEUNIT')
    wavelnth_comment = header['KEYCOMMENTS'].get('WAVELNTH')
    waveunit = header.get('WAVEUNIT')
    if waveunit is not None:
        metre_submultiples = {
            0: parse_waveunit_comment(waveunit_comment),
            -1: 'dm',
            -2: 'cm',
            -3: 'mm',
            -6: 'um',
            -9: 'nm',
            -10: 'angstrom',
            -12: 'pm',
            -15: 'fm',
            -18: 'am',
            -21: 'zm',
            -24: 'ym'}
        waveunit = metre_submultiples.get(waveunit, str(waveunit).lower())
    elif waveunit_comment is not None:
        waveunit = parse_waveunit_comment(waveunit_comment)
    elif wavelnth_comment is not None:
        # supported formats (where $UNIT is the unit like "nm" or "Angstrom"):
        #   "Observed wavelength ($UNIT)"
        #   "[$UNIT] ..."
        parentheses_pattern = r'Observed wavelength \((\w+?)\)$'
        brackets_pattern = r'^\[(\w+?)\]'
        for pattern in [parentheses_pattern, brackets_pattern]:
            m = re.search(pattern, wavelnth_comment)
            if m is not None:
                waveunit = m.group(1)
                break
    if waveunit == '':
        return None # To fix problems associated with HMI FITS.
    return waveunit
