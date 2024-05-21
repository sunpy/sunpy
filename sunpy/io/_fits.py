"""
This module provides a FITS file reader for internal use.

Notes
-----
1. FITS files allow comments to be attached to every value in the header.
   This is implemented in this module as a KEYCOMMENTS dictionary in the
   sunpy header. To add a comment to the file on write, add a comment to this
   dictionary with the same name as a key in the header (upcased).

2. Due to the way `~astropy.io.fits` works with images, the header dictionary may
   differ depending on whether is accessed before or after the fits[0].data
   is requested. If the header is read before the data then the original
   header will be returned. If the header is read after the data has been
   accessed then the data will have been scaled and a modified header
   reflecting these changes will be returned: BITPIX may differ and
   BSCALE and B_ZERO may be dropped in the modified version.

3. The verify('silentfix+warn') call attempts to handle violations of the FITS
   standard. For example, ``nan`` values will be converted to "nan" strings.
   Attempting to cast a `astropy.io.fits.Header` to a dictionary while it contains
   invalid header tags will result in an error so verifying it early on
   makes the header easier to work with later.
"""
import os
import re
import sys
import math
import traceback
import collections
import collections.abc

from astropy.io import fits

from sunpy.io._header import FileHeader
from sunpy.util.exceptions import warn_metadata, warn_user
from sunpy.util.io import HDPair

__all__ = ['header_to_fits', 'read', 'get_header', 'write', 'extract_waveunit', 'format_comments_and_history']


def read(filepath, hdus=None, memmap=None, **kwargs):
    """
    Read a fits file.

    Parameters
    ----------
    filepath : `str`
        The fits file to be read.
    hdus : `int` or iterable
        The HDU indexes to read from the file.
    **kwargs : `dict`, optional
        Passed to `astropy.io.fits.open`.

    Returns
    -------
    `list`
        A list of (data, header) tuples

    Notes
    -----
    This routine reads all the HDU's in a fits file and returns a list of the
    data and a FileHeader instance for each one.

    Also all comments in the original file are concatenated into a single
    "comment" key in the returned FileHeader.
    """
    with fits.open(filepath, ignore_blank=True, memmap=memmap, **kwargs) as hdulist:
        if hdus is not None:
            if isinstance(hdus, int):
                hdulist = hdulist[hdus]
            elif isinstance(hdus, collections.abc.Iterable):
                hdulist = [hdulist[i] for i in hdus]

        hdulist = fits.hdu.HDUList(hdulist)
        for h in hdulist:
            h.verify('silentfix+warn')

        headers = get_header(hdulist)
        pairs = []

        for i, (hdu, header) in enumerate(zip(hdulist, headers)):
            try:
                pairs.append(HDPair(hdu.data, header))
            except (KeyError, ValueError) as e:
                message = f"Error when reading HDU {i}. Skipping.\n"
                for line in traceback.format_tb(sys.exc_info()[2]):
                    message += line
                    message += '\n'
                message += repr(e)
                warn_user(message)

    return pairs


def get_header(afile):
    """
    Read a fits file and return just the headers for all HDU's.

    Parameters
    ----------
    afile : `str` or `astropy.io.fits.HDUList`
        The file to be read, or HDUList to process.

    Returns
    -------
    `list`
        A list of `sunpy.io._header.FileHeader` headers.
    """
    if isinstance(afile, fits.HDUList):
        hdulist = afile
        close = False
    else:
        hdulist = fits.open(afile, ignore_blank=True)
        hdulist.verify('silentfix')
        close = True

    try:
        headers = []
        for hdu in hdulist:
            headers.append(format_comments_and_history(hdu.header))
    finally:
        if close:
            hdulist.close()
    return headers


def format_comments_and_history(input_header):
    """
    Combine ``COMMENT`` and ``HISTORY`` cards into single
    entries. Extract ``KEYCOMMENTS`` into a single entry
    and put ``WAVEUNIT`` into its own entry.

    Parameters
    ----------
    input_header : `~astropy.io.fits.Header`
        The header to be processed.

    Returns
    -------
    `sunpy.io._header.FileHeader`
    """
    try:
        comment = "\n".join(input_header['COMMENT']).strip()
    except KeyError:
        comment = ""
    try:
        history = "\n".join(input_header['HISTORY']).strip()
    except KeyError:
        history = ""

    header = FileHeader(input_header)
    header['COMMENT'] = comment
    header['HISTORY'] = history

    # Strip out KEYCOMMENTS to a dict, the hard way
    keydict = {}
    for card in input_header.cards:
        if card.comment != '':
            keydict.update({card.keyword: card.comment})
    header['KEYCOMMENTS'] = keydict
    waveunit = extract_waveunit(header)
    if waveunit is not None:
        header['WAVEUNIT'] = waveunit
    return header


def write(fname, data, header, hdu_type=None, **kwargs):
    """
    Take a data header pair and write a FITS file.

    Parameters
    ----------
    fname : `str`
        File name, with extension.
    data : `numpy.ndarray`
        n-dimensional data array.
    header : `dict`
        A header dictionary.
    hdu_type : `~astropy.io.fits.hdu.base.ExtensionHDU` instance or class, optional
        By default, a FITS file is written with the map in its primary HDU.
        If a type is given, a new HDU of this type will be created.
        If a HDU instance is given, its data and header will be updated from the map.
        Then that HDU instance will be written to the file.
    kwargs :
        Additional keyword arguments are given to
        `~astropy.io.fits.HDUList.writeto`.
    """
    # Copy header so the one in memory is left alone
    header = header.copy()

    fits_header = header_to_fits(header)

    if isinstance(fname, str):
        fname = os.path.expanduser(fname)

    fitskwargs = {'output_verify': 'fix'}
    fitskwargs.update(kwargs)

    if not hdu_type:
        hdu_type = fits.PrimaryHDU

    if isinstance(hdu_type, fits.PrimaryHDU | fits.hdu.base.ExtensionHDU):
        hdu = hdu_type  # HDU already initialized
        # Merge `header` into HDU's header
        # Values in `header` take priority, including cards such as
        # 'SIMPLE' and 'BITPIX'.
        hdu.header.extend(fits_header, strip=False, update=True)
        # Set the HDU's data
        hdu.data = data
    else:
        hdu = hdu_type(data=data, header=fits_header)

    if not isinstance(hdu, fits.PrimaryHDU):
        hdul = fits.HDUList([fits.PrimaryHDU(), hdu])
    else:
        hdul = fits.HDUList([hdu])

    hdul.writeto(fname, **fitskwargs)


def header_to_fits(header):
    """
    Convert a header dict to a `~astropy.io.fits.Header`.
    """
    # Copy the header to avoid modifying it in place
    header = header.copy()
    # The comments need to be added to the header separately from the normal
    # kwargs. Find and deal with them:
    fits_header = fits.Header()
    # Check Header
    key_comments = header.pop('KEYCOMMENTS', False)
    # This first iteration separates out COMMENTS and HISTORY into their own
    # entries.
    header_items = []
    for k, v in header.items():
        if k.upper() in ('COMMENT', 'HV_COMMENT', 'HISTORY'):
            for v_line in str(v).split('\n'):
                header_items.append((k, v_line))
        else:
            header_items.append((k, v))

    for k, v in header_items:
        # Drop any keys that have non-ascii characters
        if not fits.Card._ascii_text_re.match(str(v)):
            warn_metadata(f'The meta key {k} is not valid ascii, dropping from the FITS header')
            continue
        # Drop any keys which are too long to save into FITS
        if len(k) > 8:
            warn_metadata(f"The meta key {k} is too long, dropping from the FITS header "
                          "(maximum allowed key length is 8 characters).")
            continue

        if isinstance(v, float) and math.isnan(v):
            warn_metadata(f'The meta key {k} has a NaN value, which is not valid in a FITS '
                          'header, dropping from the FITS header')
            continue

        if k.upper() in ('COMMENT', 'HV_COMMENT'):
            fits_header.add_comment(v)
        elif k.upper() == 'HISTORY':
            fits_header.add_history(v)
        elif isinstance(v, fits.header._HeaderCommentaryCards):
            if k != '':
                fits_header.append(fits.Card(k, str(v).split('\n')))
        else:
            # For some horrific reason, we save a list to the wavelnth key in
            # sources/rhessi.py. This is the least invasive fix for that stupidity.
            if isinstance(v, list):
                v = str(v)
            fits_header.append(fits.Card(k, v))

    if isinstance(key_comments, dict):
        for k, v in key_comments.items():
            # Check that the Card for the comment exists before trying to write to it.
            if k in fits_header:
                fits_header.comments[k] = v
    elif key_comments:
        raise TypeError("KEYCOMMENTS must be a dictionary")

    return fits_header


def extract_waveunit(header):
    """
    Attempt to read the wavelength unit from a given FITS header.

    Parameters
    ----------
    header : `sunpy.io._header.FileHeader`
        One `~sunpy.io._header.FileHeader` instance which was created by
        reading a FITS file. For example, `sunpy.io._fits.get_header` returns a list of
        such instances.

    Returns
    -------
    `str`
        The wavelength unit that could be found or ``None`` otherwise.

    Examples
    --------
    The goal of this function is to return a string that can be used in
    conjunction with the astropy.units module so that the return value can be
    directly passed to `astropy.units.Unit`.

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

    if 'KEYCOMMENTS' not in header:
        return None
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
        #   "Wavelength of obs ($UNIT)"
        #   "[$UNIT] ..."
        patterns = [
            r'Observed wavelength \((\w+?)\)$',
            r'Wavelength of obs \((\w+?)\)$',
            r'^\[(\w+?)\]',
        ]
        for pattern in patterns:
            m = re.search(pattern, wavelnth_comment)
            if m is not None:
                waveunit = m.group(1)
                break
    if waveunit == '':
        return None  # To fix problems associated with HMI FITS.
    return waveunit
