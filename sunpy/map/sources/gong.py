"""
GONG Map subclass definitions
"""
import numpy as np
from astropy.time import Time

from sunpy.map import GenericMap

__all__ = ['GongSynopticMap']


class GongSynopticMap(GenericMap):
    def __init__(self, data, header, **kwargs):
        # Fix coordinate system stuff
        if 'KEYCOMMENTS' in header:
            if 'deg' in header['KEYCOMMENTS']['CDELT1']:
                header['CUNIT1'] = 'deg'
            if header['KEYCOMMENTS']['CDELT2'] == 'Sine-lat step':
                header['CUNIT2'] = 'deg'
                # Instead of the spacing in sin(lat), this should be 180/pi times
                # that value (see Thompson 2005)
                header['CDELT2'] = 180 / np.pi * header['CDELT2']
                header['KEYCOMMENTS']['CDELT2'] = '180 / pi * sine-lat step'
        # Fix timestamp
        if 'time-obs' in header:
            header['date-obs'] = (header['date-obs'] + ' ' +
                                  header.pop('time-obs'))
            header['date-obs'] = Time(header['date-obs']).isot
        # Fix unit
        if 'bunit' in header and header['bunit'] == 'Gauss':
            header['bunit'] = 'G'
        # Fix observer coordinate
        if 'hglt_obs' not in header:
            header.update(_earth_obs_coord_meta(header['date-obs']))
        super().__init__(data, header, **kwargs)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an GONG map."""
        return (str(header.get('TELESCOP', '')).endswith('GONG') and
                str(header.get('CTYPE1', '').startswith('CRLN')))
