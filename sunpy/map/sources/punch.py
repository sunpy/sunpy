"""PUNCH Map subclass definitions"""

from matplotlib.colors import LogNorm

import astropy.units as u

from sunpy.map.mapbase import GenericMap

__all__ = ['PUNCHMap']

class PUNCHMap(GenericMap):
    """
    PUNCH Map.

    The Polarimeter to Unify the Corona and Heliosphere (PUNCH) is a constellation observatory consisting of four satellites, together imaging the solar corona in polarized white light.
    Coordinated observations from three Wide Field Imagers (WFIs) and one Near Field Imager (NFI) are processed and meshed into one virtual observatory with a field of view spanning 45-degrees in radius away from the Sun.

    PUNCH launched and began operations on 11 March 2025.

    Notes
    -----

    References
    ----------
    * `PUNCH Mission Page <https://punch.space.swri.edu>`__
    * `PUNCH SOC Development <https://github.com/punch-mission/>`__
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        self.nickname = f"{self.observatory} - {self.instrument}"
        self.plot_settings["cmap"] = "punch"

        self.plot_settings["norm"] = LogNorm()


    @property
    def unit(self):
        """Units of the observation."""
        unit_str = self.meta.get('bunit', self.meta.get("BUNIT"))
        if unit_str is None:
            return
        return u.def_unit(unit_str)

    # Used by the Map factory to determine if this subclass should be used
    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if data, header corresponds to a PUNCH image."""
        # Returns True only if this is data and header from PUNCH
        instrument = str(header.get('instrume', ''))
        return instrument.startswith('NFI') or instrument.startswith('WFI') or header.get('obsrvtry', '').startswith('PUNCH')
