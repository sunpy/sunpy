"""RHESSI Map subclass definitions"""

__author__ = "Steven Christe"
__email__ = "steven.d.christe@nasa.gov"

from sunpy.map import GenericMap

__all__ = ['RHESSIMap']


class RHESSIMap(GenericMap):
    """RHESSI Image Map.

    The RHESSI mission consists of a single spin-stabilized
    spacecraft in a low-altitude orbit inclined 38 degrees to
    the Earth's equator. The only instrument on board is an
    Germaniun imaging spectrometer with the ability to obtain high
    fidelity solar images in X rays (down to 3 keV) to gamma rays (1 MeV).

    RHESSI provides an angular resolution of 2 arcseconds at
    X-ray energies below ~40 keV, 7 arcseconds to 400 keV,
    and 36 arcseconds for gamma-ray lines and continuum above 1 MeV.

    RHESSI was launched on 5 February 2002.

    References
    ----------
    * RHESSI Homepage `<https://hesperia.gsfc.nasa.gov/rhessi3/index.html>`_
    * Mission Paper `<https://doi.org/10.1023/A:1022428818870>`_

    .. warning::

        This software is in beta and cannot read fits files containing more than one image.
    """

    def __init__(self, data, header, **kwargs):
        # Fix some broken/misapplied keywords
        if header['ctype1'] == 'arcsec':
            header['cunit1'] = 'arcsec'
            self._fix_and_warn_header(header, 'ctype1', 'HPLN-TAN', replace_old=True)
        if header['ctype2'] == 'arcsec':
            header['cunit2'] = 'arcsec'
            self._fix_and_warn_header(header, 'ctype2', 'HPLT-TAN', replace_old=True)
        super().__init__(data, header, **kwargs)

        self._nickname = self.detector
        # TODO Currently (8/29/2011), cannot read fits files containing more
        # than one image (schriste)
        self._fix_and_warn('waveunit', 'keV')
        self.meta['wavelnth'] = [self.meta['energy_l'], self.meta['energy_h']]
        self.plot_settings['cmap'] = 'rhessi'

    @property
    def detector(self):
        """
        Returns the name of the detector
        """
        return self.meta['telescop']

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an RHESSI image"""
        return header.get('instrume') == 'RHESSI'
