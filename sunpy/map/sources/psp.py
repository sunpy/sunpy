"""
Parker Solar Probe subclass definitions.
"""

from sunpy.map import GenericMap

__all__ = ['WISPRMap']


class WISPRMap(GenericMap):
    """
    WISPR Map

    The The Wide-field Imager for Parker Solar Probe (WISPR) is a white light
    telescope onboard the Parker Solar Probe (PSP) spacecraft.

    References
    ----------
    * `PSP science gateway <https://sppgway.jhuapl.edu//>`__
    * `WISPR Instrument Page <https://wispr.nrl.navy.mil//>`__
    * `Instrument Paper <https://doi.org/10.1007/s11214-014-0114-y>`__
    """
    @property
    def processing_level(self):
        lvl = self.meta.get('level', None)
        if lvl is None:
            return
        return int(lvl[1])

    @property
    def exposure_time(self):
        return self.meta.get('xposure', 0.0) * self.timeunit

    @property
    def name(self):
        return 'WISPR ' + super().name

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an WISPR image"""
        is_psp = 'parker solar probe' in str(header.get('obsrvtry', '')).lower()
        is_wispr = str(header.get('instrume', '')).startswith('WISPR')
        return is_psp and is_wispr
