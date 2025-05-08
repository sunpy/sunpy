"""Yohkoh SXT Map subclass definitions"""

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"


import astropy.units as u
from astropy.coordinates import ITRS, SphericalRepresentation
from astropy.visualization import PowerStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch

__all__ = ['SXTMap']


class SXTMap(GenericMap):
    """Yohkoh SXT Image Map

    The Yohkoh Soft X-ray Telescope (SXT) observed the full solar disk
    (42 x 42 arcminutes) in the 0.25 - 4.0 keV range.
    It consisted of a glancing incidence mirror and a CCD sensor and
    used thin metallic filters to acquire images in restricted
    portions of its energy range.
    SXT could resolve features down to 2.5 arcseconds.
    Information about the temperature and density of the plasma
    emitting the observed x-rays was obtained by comparing images acquired with
    the different filters.
    Images could be obtained every 2 to 8 seconds.
    Smaller images with a single filter could be obtained as frequently as
    once every 0.5 seconds.
    Yohkoh was launched on 30 August 1991 and ceased operations on
    14 December 2001.

    Notes
    -----
    Observer location:  We use the ITRS coordinates provided in
    the FITS header (``LAT``, ``LONG``, ``RADIUS``)
    for the spacecraft location even when coordinates in the heliographic
    Stonyhurst frame are provided (``HGS_LON``, ``HGS_LAT``, ``DSUN_OBS``).
    The two locations differ substantially, and the ITRS coordinates
    are more accurate than the HGS information.

    Using this information for the observer coordinate was implemented
    in ``sunpy`` 4.0.6 to work with the FITS files that were re-processed in 2016.

    References
    ----------
    * `Yohkoh Mission Page <http://solar.physics.montana.edu/sxt/>`__
    * `Fits header reference <http://proba2.oma.be/data/SWAP/level0>`__
    * `Yohkoh Analysis Guide <http://ylstone.physics.montana.edu/ylegacy/yag.html>`__
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)

        self.plot_settings['cmap'] = 'yohkohsxt' + self.measurement[0:2].lower()
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, PowerStretch(0.5)), clip=False)

    @property
    def _supported_observer_coordinates(self):
        return [(('long', 'lat', 'radius'), {'lon': self.meta.get('long'),
                                             'lat': self.meta.get('lat'),
                                             'distance': self.meta.get('radius'),
                                             'unit': (u.deg, u.deg, u.km),
                                             'representation_type': SphericalRepresentation,
                                             'frame': ITRS, })
                ] + super()._supported_observer_coordinates

    @property
    def observatory(self):
        return "Yohkoh"

    @property
    def detector(self):
        return "SXT"

    @property
    def measurement(self):
        s = self.meta.get('wavelnth', '')
        if s == 'Al.1':
            s = 'Al01'
        elif s.lower() == 'open':
            s = 'white-light'
        return s

    @property
    def wavelength(self):
        """
        Returns `None`, as SXT is a broadband imager.
        """

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an SXT image"""
        return header.get('instrume') == 'SXT'
