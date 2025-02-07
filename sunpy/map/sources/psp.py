"""
Parker Solar Probe subclass definitions.
"""

from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.map import GenericMap

__all__ = ['WISPRMap']


class WISPRMap(GenericMap):
    """
    WISPR Map

    The The Wide-field Imager for Parker Solar Probe (WISPR) is a white light
    telescope onboard the Parker Solar Probe (PSP) spacecraft.

    Notes
    -----
    By default, plotting of this map will set the lower bound to zero
    (i.e., clip out negative values for pixels). You can change this bound
    by modifying ``.plot_settings['norm'].vmin``.

    References
    ----------
    * `PSP science gateway <https://sppgway.jhuapl.edu//>`__
    * `WISPR Instrument Page <https://wispr.nrl.navy.mil//>`__
    * :cite:t:`vourlidas_wide-field_2016`
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)

        self.plot_settings['norm'] = ImageNormalize(
            stretch=AsinhStretch(a=0.001), vmin=0)

    @property
    def processing_level(self):
        lvl = self.meta.get('level', None)
        if lvl is None:
            return
        # Chop off the leading 'L' if present
        if lvl[0] == 'L':
            lvl = lvl[1:]
        try:
            lvl = int(lvl)
        except ValueError:
            # The int conversion will fail for L2b files, and we should fail
            # safe if the user chooses to customize this with other values.
            pass
        return lvl

    @property
    def name(self):
        return 'WISPR ' + super().name

    @property
    def detector(self):
        detector = self.meta.get('detector', "")
        if detector == 1:
            return "Inner"
        if detector == 2:
            return "Outer"
        # Official data products will only be 1 or 2, but we should fail safe
        # if users customize this value themselves.
        return detector

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an WISPR image"""
        is_psp = 'parker solar probe' in str(header.get('obsrvtry', '')).lower()
        is_wispr = str(header.get('instrume', '')).startswith('WISPR')
        return is_psp and is_wispr
