import astropy.units as u

from sunpy.map.mapbase import GenericMap, SpatialPair

__all__ = ['RHESSIMap']


class RHESSIMap(GenericMap):
    """
    RHESSI Image Map.

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
    * RHESSI Homepage `<https://hesperia.gsfc.nasa.gov/rhessi3/index.html>`__
    * :cite:t:`lin_reuven_2002`

    .. warning::

        Cannot read fits files containing more than one image.
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        self._nickname = self.detector
        self.plot_settings['cmap'] = 'rhessi'

    def _get_cmap_name(self):
        return "rhessi"

    @property
    def _timesys(self):
        """
        RHESSI maps can incorrectly use the TIMESYS keyword for the reference
        time. If this is the case, returns the FITS default UTC.
        """
        if ('TIMESYS' in self.meta and
                self.meta['keycomments']['TIMESYS'] == 'Reference Time'):
            return 'UTC'
        else:
            return super()._timesys

    def _rotation_matrix_from_crota(self):
        """
        RHESSI maps can have their rotation in CROTA.
        """
        return super()._rotation_matrix_from_crota(crota_key='CROTA')

    @property
    def spatial_units(self):
        """
        If CTYPE{1 or 2} are equal to 'arcsec', assumes that CTYPE{1 or 2}
        respectively are intended to be 'arcsec'.
        """
        units = [self.meta.get('cunit1', None), self.meta.get('cunit2', None)]
        if self.meta['ctype1'] == 'arcsec':
            units[0] = 'arcsec'
        if self.meta['ctype2'] == 'arcsec':
            units[1] = 'arcsec'
        units = [None if unit is None else u.Unit(unit.lower()) for unit in units]
        return SpatialPair(units[0], units[1])

    @property
    def coordinate_system(self):
        """
        If CTYPE{1 or 2} are equal to 'arcsec', assumes that CTYPE{1 or 2}
        respectively are intended to be 'HPLN-TAN' or 'HPLT-TAN'.
        """
        ctype1, ctype2 = self.meta['ctype1'], self.meta['ctype2']
        if ctype1 == 'arcsec':
            ctype1 = 'HPLN-TAN'
        if ctype2 == 'arcsec':
            ctype2 = 'HPLT-TAN'
        return SpatialPair(ctype1, ctype2)

    @property
    def waveunit(self):
        """
        If the WAVEUNIT FITS keyword is not present, defaults to keV.
        """
        unit = self.meta.get("waveunit", 'keV')
        return u.Unit(unit)

    @property
    def wavelength(self):
        return u.Quantity([self.meta['energy_l'], self.meta['energy_h']],
                          unit=self.waveunit)

    @property
    def detector(self):
        return self.meta['telescop']

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an RHESSI image"""
        return header.get('instrume') == 'RHESSI'
