import astropy.units as u

from sunpy.map import GenericMap

__all__ = ['SJIMap']


class SJIMap(GenericMap):
    """
    A 2D IRIS Slit Jaw Imager Map.

    The Interface Region Imaging Spectrograph (IRIS) small explorer spacecraft
    provides simultaneous spectra and images of the photosphere, chromosphere,
    transition region, and corona with 0.33 to 0.4 arcsec spatial resolution,
    2-second temporal resolution and 1 km/s velocity resolution over a
    field-of- view of up to 175 arcsec by 175 arcsec.  IRIS consists of a 19-cm
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
    * `IRIS Mission Page <https://iris.lmsal.com>`_
    * `IRIS Analysis Guide <https://iris.lmsal.com/itn26/itn26.pdf>`_
    * `IRIS Instrument Paper <https://doi.org/10.1007/s11207-014-0485-y>`_
    """

    def __init__(self, data, header, **kwargs):
        # Assume pixel units are arcesc if not given
        header['cunit1'] = header.get('cunit1', 'arcsec')
        header['cunit2'] = header.get('cunit2', 'arcsec')
        header['bunit'] = header.get('bunit', 'counts').lower().replace('corrected', '').strip().upper()
        super().__init__(data, header, **kwargs)

        self.meta['detector'] = "SJI"
        self.meta['waveunit'] = u.Unit(self.meta.get('waveunit', "Angstrom"))
        self.meta['wavelnth'] = self.meta.get('wavelnth', self.meta.get('twave1')) * self.waveunit

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an IRIS SJI image"""
        tele = str(header.get('TELESCOP', '')).startswith('IRIS')
        obs = str(header.get('INSTRUME', '')).startswith('SJI')
        return tele and obs
