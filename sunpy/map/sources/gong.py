"""
GONG Map subclass definitions
"""
import numpy as np
from astropy.time import Time
import astropy.units as u

from sunpy.map import GenericMap
from sunpy.coordinates import get_earth, HeliographicStonyhurst

__all__ = ['GongSynopticMap', 'ADAPTMap']


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


class ADAPTMap(GenericMap):
    def __init__(self, data, header, **kwargs):
        if 'date-obs' not in header:
            header['date-obs'] = header['maptime']
        # Fix CTYPE
        if header['ctype1'] == 'Long':
            header['ctype1'] = 'CRLN-CAR'
        if header['ctype2'] == 'Lat':
            header['ctype2'] = 'CRLT-CAR'

        super().__init__(data, header, **kwargs)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an ADAPT map."""
        return header.get('model') == 'ADAPT'


def _observer_coord_meta(observer_coord):
    """
    Convert an observer coordinate into FITS metadata.
    """
    new_obs_frame = HeliographicStonyhurst(
        obstime=observer_coord.obstime)
    observer_coord = observer_coord.transform_to(new_obs_frame)

    new_meta = {}
    new_meta['hglt_obs'] = observer_coord.lat.to_value(u.deg)
    new_meta['hgln_obs'] = observer_coord.lon.to_value(u.deg)
    new_meta['dsun_obs'] = observer_coord.radius.to_value(u.m)
    return new_meta


def _earth_obs_coord_meta(obstime):
    """
    Return metadata for an Earth obeserver coordinate.
    """
    return _observer_coord_meta(get_earth(obstime))
