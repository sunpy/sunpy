# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
'''EIS spectral cube definitions'''

from __future__ import absolute_import

from astropy.io import fits
from astropy import wcs
from sunpy.cube import SpectralCube
import re

__all__ = ['EISSpectralCube']


def _clean(header):
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


class EISSpectralCube(SpectralCube):
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
        wcs: astropy.wcs.WCS object
            The world coordinate system for the array.
        window: int
            The window this cube belongs to in the file. Used to fetch the
            correct metadata from the header
        dataHeader: astropy.io.fits.Header object
            The header for the BINTableHDU section of the FITS file
        primaryHeader: astropy.io.fits.Header object
            The main header for the whole file.
        '''
        h = _dictionarize_header(dataHeader, primaryHeader, window)
        SpectralCube.__init__(self, data.T, wcs, meta=h)
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
        # TODO: Make sure each cube ahas a correct wcs.
        w = wcs.WCS(header=header, naxis=3)
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


def _isInWindow(key, window):
    l = re.findall(r'\d+$', key)  # finds numbers at the end of the key
    if len(l) == 0:
        return False
    else:
        return window == int(l[0])


def _dictionarize_header(dataHeader, primaryHeader, window):
    ph = dict(primaryHeader)
    dh = {}
    for k in dataHeader:
        if _isInWindow(k, window):
            newkey = re.sub(r'\d+$', '', k)
            dh[newkey] = dataHeader[k]

    ph.update(dh)
    return ph
