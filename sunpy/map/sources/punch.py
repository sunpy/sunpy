"""PUNCH Map subclass definitions"""

from astropy.visualization import ImageNormalize, LogStretch

from sunpy.map.mapbase import GenericMap
from sunpy.map.sources.source_type import source_stretch

__all__ = ['PUNCHMap']

class PUNCHMap(GenericMap):
    """
    PUNCH Map.

    The Polarimeter to Unify the Corona and Heliosphere (PUNCH) is a constellation observatory
    consisting of four satellites, together imaging the solar corona in polarized white light.
    Coordinated observations from three Wide Field Imagers (WFIs) and one Near Field Imager (NFI)
    are processed and meshed into one virtual observatory with a field of view spanning
    45-degrees in radius away from the Sun.

    PUNCH launched and began operations on 11 March 2025.

    References
    ----------
    * `PUNCH Mission Page <https://punch.space.swri.edu>`__
    * `PUNCH SOC Development <https://github.com/punch-mission/>`__
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        self.nickname = f"{self.observatory} - {self.instrument}"
        self.plot_settings["cmap"] = "punch"
        self.plot_settings["norm"] = ImageNormalize(
            stretch=source_stretch(self.meta, LogStretch()), clip=False)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if data, header corresponds to a PUNCH image."""
        return header.get('obsrvtry', '').startswith('PUNCH')
