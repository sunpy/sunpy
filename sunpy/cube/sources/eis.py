# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
'''EIS spectral cube definitions'''

from __future__ import absolute_import

from astropy.io import fits
from sunpy.wcs.wcs import WCS
from sunpy.cube import Cube
import re

__all__ = ['EISSpectralCube']


def _clean(header):
    # TODO: find a way to identify cubes containing time
    """ Fixes non-standard or deprecated CTYPEn FITS keywords.

    Parameters
    ----------
    header : astropy.io.fits.Header
        The header to be cleaned.
    """
    header['ctype1'] = 'HPLN-TAN'  # Helioprojective longitude, TAN projection
    header['ctype2'] = 'HPLT-TAN'  # Helioprojective latitude, TAN projection
    header['ctype3'] = 'WAVE   '  # Wavelength axis, default (TAB) projection
    header['naxis'] = 3
    return header


class EISSpectralCube(Cube):
    '''EIS Spectral Cube subclass.

    References
    ----------
    For an overview of the mission
    http://solarb.mssl.ucl.ac.uk/SolarB/
    '''
    def __init__(self,
                 data, wcs, window=1, dataHeader=None, primaryHeader=None):
        '''
        Constructor function.

        Parameters
        ----------
        data: numpy ndarray
            The cube containing the data
        wcs: sunpy.wcs.wcs.WCS object
            The world coordinate system for the array.
        window: int
            The window this cube belongs to in the file. Used to fetch the
            correct metadata from the header
        dataHeader: astropy.io.fits.Header object
            The header for the BINTableHDU section of the FITS file
        primaryHeader: astropy.io.fits.Header object
            The main header for the whole file.
        '''
        header = _dictionarize_header(dataHeader, primaryHeader, window)
        Cube.__init__(self, data.T, wcs, meta=header)
        # Data is transposed here because EIS orders (y, lambda) by x or time,
        # not (y, x) by lambda.

    @classmethod
    def read(cls, filename, **kwargs):
        """ Reads in a given FITS file and returns a dictionary of new
        EISSpectralCubes. Additional parameters are given to fits.open.

        Parameters
        ----------
        filename : string
            Complete location of the FITS file
        """
        hdulist = fits.open(name=filename, **kwargs)
        header = _clean(hdulist[0].header)
        # TODO: Make sure each cube has a correct wcs.
        w = WCS(header=header, naxis=3)
        wavelengths = [c.name for c in hdulist[1].columns if c.dim is not None]
        data = [hdulist[1].data[wav] for wav in wavelengths]
        cubes = [EISSpectralCube(data[i], w, i+1, dataHeader=hdulist[1].header,
                                 primaryHeader=hdulist[0].header)
                 for i in range(len(data))]
        return dict(zip(wavelengths, cubes))

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        # TODO: This will have to be changed once other sources are added.
        return True


def _is_in_window(key, window):
    '''
    Checks if a given key forms part of the specified spectral window.

    Parameters
    ----------
    key: str
        The key to be validated
    window: int
        The desired window
    '''
    end = re.findall(r'\d+$', key)  # finds numbers at the end of the key
    if len(end) == 0:
        return False
    else:
        return window == int(end[0])


def _dictionarize_header(data_header, primary_header, window):
    '''
    Combines the given FITS primary header and the bintable header for a
    specified window into a dictionary.

    Parameters
    ----------
    data_header: astropy.io.fits.Header object, dict, or dict-like object.
        secondary header to be pruned for the specified window
    primary_header: astropy.io.fits.Header object, dict, or dict-like object.
        The main FITS file header
    window: int
        The window to be chosen out of the data header.
    '''
    ph = dict(primary_header)
    dh = {}
    for k in data_header:
        if _is_in_window(k, window):
            newkey = re.sub(r'\d+$', '', k)
            dh[newkey] = data_header[k]

    ph.update(dh)
    return ph
