
import astropy.units as u

from sunpy.map.mapbase import GenericMap, SpatialPair

__all__ = ['SJIMap']


class SJIMap(GenericMap):
    """
    A 2D IRIS Slit Jaw Imager Map.

    The Interface Region Imaging Spectrograph (IRIS) small explorer spacecraft
    provides simultaneous spectra and images of the photosphere, chromosphere,
    transition region, and corona with 0.33 to 0.4 arcsec spatial resolution,
    2-second temporal resolution and 1 km/s velocity resolution over a
    field-of- view of up to 175 arcsec by 175 arcsec. IRIS consists of a 19-cm
    UV telescope that feeds a slit-based dual-bandpass imaging spectrograph.

    Slit-jaw images in four different passbands (C ii 1330, Si iv 1400,
    Mg ii k 2796 and Mg ii wing 2830  A) can be taken simultaneously with
    spectral rasters that sample regions up to 130 arcsec by 175 arcsec at a
    variety of spatial samplings (from 0.33 arcsec and up).
    IRIS is sensitive to emission from plasma at temperatures between
    5000 K and 10 MK.

    IRIS was launched into a Sun-synchronous orbit on 27 June 2013.

    References
    ----------
    * `IRIS Mission Page <https://iris.lmsal.com>`__
    * `IRIS Analysis Guide <https://iris.lmsal.com/itn26/itn26.pdf>`__
    * :cite:t:`de_pontieu_interface_2014`
    """
    @property
    def detector(self):
        return "SJI"

    @property
    def spatial_units(self):
        """
        If not present in CUNIT{1,2} keywords, defaults to arcsec.
        """
        return SpatialPair(u.Unit(self.meta.get('cunit1', 'arcsec')),
                           u.Unit(self.meta.get('cunit2', 'arcsec')))

    @property
    def waveunit(self):
        """
        Taken from WAVEUNIT, or if not present defaults to Angstrom.
        """
        return u.Unit(self.meta.get('waveunit', "Angstrom"))

    @property
    def wavelength(self):
        """
        Taken from WAVELNTH, or if not present TWAVE1.
        """
        return self.meta.get('wavelnth', self.meta.get('twave1')) * self.waveunit

    @property
    def unit(self):
        unit_str = self.meta.get('bunit', None)
        if unit_str is None:
            return
        # Remove "corrected" so that the unit can be parsed
        unit_str = unit_str.lower().replace('corrected', '').strip()
        return self._parse_fits_unit(unit_str)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an IRIS SJI image"""
        tele = str(header.get('TELESCOP', '')).startswith('IRIS')
        obs = str(header.get('INSTRUME', '')).startswith('SJI')
        return tele and obs
