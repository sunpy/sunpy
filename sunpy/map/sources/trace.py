"""TRACE Map subclass definitions"""
from __future__ import absolute_import, division, absolute_import
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"

from sunpy.map import GenericMap
from sunpy.cm import cm

from matplotlib import colors

__all__ = ['TRACEMap']

class TRACEMap(GenericMap):
    """TRACE Image Map

    The Transition Region and Coronal Explorer was a
    NASA Small Explorer (SMEX) mission to image the
    solar corona and transition region at high angular and temporal resolution.
    TRACE observed the Sun in the following passbands, 5000 A, 1700 A, 1600 A,
    1550 A (C IV), 1216 A (H1 Lyman-alpha), 173 A (Fe IX), 195 A (Fe XII),
    and 284 A (Fe XV). TRACE provides solar images with an 8.5 x 8.5 arcminute
    field of view and 0.5 arcsecond pixels. It was placed in a sun-synchronous
    orbit, enabling it to make continuous solar observations.

    The TRACE mission operated was launched on 2 April 1998 and obtained its
    last science image on 6 June 2010 23:56 UT.

    References
    ----------
    * `Mission/Instrument Page <http://trace.lmsal.com>`_
    * `Fits headers <http://trace.lmsal.com/Project/Instrument/cal/>`_
    * `Analysis Guide <http://trace.lmsal.com/tag/>`_
    * `Passband reference <http://trace.lmsal.com/Project/Instrument/inspass.htm>`_

    .. note::

        Note that this map definition is currently only being tested on JPEG2000
        files. TRACE FITS data is stored in a more complex format. Typically
        TRACE data is stored in hourly "tri" files that store all the data taken
        by TRACE in the hour indicated by the filename. Those files must first be
        understood and parsed to obtain the science data. The ability to do this
        is not yet in SunPy, but is available in SSWIDL. Please refer to the links
        above concerning how to read "tri" files in SSWIDL.
    """

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        # It needs to be verified that these must actually be set and are not
        # already in the header.
        self.meta['detector'] = "TRACE"
        self.meta['obsrvtry'] = "TRACE"
        self._nickname = self.detector
        # Colour maps
        self.plot_settings['cmap'] = cm.get_cmap('trace' + self.measurement)
        self.plot_settings['norm'] = colors.LogNorm()

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an TRACE image"""
        return header.get('instrume') == 'TRACE'
