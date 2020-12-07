"""Solar Orbiter Map subclass definitions"""
from astropy.visualization import AsinhStretch, ImageNormalize

from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch

__all__ = ['EUIMap']


class EUIMap(GenericMap):
    """EUI Image Map

    References
    ----------
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        self._nickname = self.detector
        self.plot_settings['cmap'] = self._get_cmap_name()
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, AsinhStretch(0.01)), clip=False)
        self.meta['exptime'] = self.meta.get('xposure', 0.0)
        # The level number is prepended by the letter L
        if self.meta.get('level'):
            self.meta['lvl_num'] = int(self.meta.get('level')[1:])
        # DN is not a FITS standard unit, so convert to counts
        if self.meta.get('bunit', None) == 'DN':
            self.meta['bunit'] = 'ct'
        if self.meta.get('bunit', None) == 'DN/s':
            self.meta['bunit'] = 'ct/s'

    @property
    def observatory(self):
        return 'SolO'

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an EUI image"""
        return str(header.get('instrume', '')).startswith('EUI')
