# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>

from __future__ import absolute_import

import spectral_cube as sc
import numpy as np
from astropy.io import fits
from astropy import wcs
from sunpy.cube import SpectralCube

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
    header['ctype3'] = 'WAVE-   '  # Wavelength axis, default (TAB) projection
    header['naxis'] = 3
    return header


class EISSpectralCube(SpectralCube):
    # TODO: write docstring
    def __init__(self, cube, dataHeader=None, primaryHeader=None):
        SpectralCube.__init__(self, cube,
                                           data_header=dataHeader,
                                           primary_header=primaryHeader)

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
        w = wcs.WCS(header=header, naxis=3)
        wavelengths = [c.name for c in hdulist[1].columns if c.dim is not None]
        data = [hdulist[1].data[wav] for wav in wavelengths]
        cubes = [sc.SpectralCube(data=d.T, wcs=w) for d in data]
        scubes = [EISSpectralCube(c, hdulist[1].header,
                                  hdulist[0].header) for c in cubes]
        return dict(zip(wavelengths, scubes))

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        return True
