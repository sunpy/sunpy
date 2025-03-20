"""PUNCH Map subclass definitions"""

import astropy.units as u

from sunpy.map.mapbase import GenericMap
from matplotlib.colors import LogNorm

__all__ = ['PUNCHMap']

class PUNCHMap(GenericMap):
    """
    PUNCH Map.

    The Polarimeter to Unify the Corona and Heliosphere (PUNCH) is a constellation observatory consisting of four satellites, together imaging the solar corona in polarized white light.

    # TODO - More mission details here - WFI, NFI, mosaic, etc.

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
        self._nickname = self.observatory + ' - ' + self.instrument
        self.plot_settings["cmap"] = "punch"
        self.plot_settings["title"] = f"{self._nickname}"
        self.plot_settings["norm"] = LogNorm(vmin=1.77e-15, vmax=3.7e-11)


    @property
    def reference_date(self):
        """
        The reference date for the coordinate system.
        """
        return self._get_date('DATE-OBS') or super().reference_date

    @property
    def unit(self):
        unit_str = self.meta.get('bunit', self.meta.get("BUNIT"))
        if unit_str is None:
            return
        return u.def_unit(unit_str)

    @property
    def observatory(self):
        return self.meta.get("OBSRVTRY")

    @property
    def instrument(self):
        return self.meta.get("INSTRUME")


    # Used by the Map factory to determine if this subclass should be used
    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """
        Determines if data, header corresponds to a PUNCH image
        """
        # Returns True only if this is data and header from PUNCH
        return header.get('instrume', '').startswith('NFI') or header.get('instrume', '').startswith('WFI') or header.get('obsrvtry', '').startswith('PUNCH')
