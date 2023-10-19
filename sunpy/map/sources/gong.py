"""
GONG Map subclass definitions
"""
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from sunpy.coordinates import HeliographicStonyhurst, get_earth
from sunpy.map import GenericMap

__all__ = ['GONGSynopticMap']

from sunpy.map.mapbase import SpatialPair


class GONGSynopticMap(GenericMap):
    """
    GONG Synoptic Map.

    The Global Oscillation Network Group (GONG) operates a six-station network of velocity
    imagers located around the Earth that observe the Sun nearly continuously. GONG
    produces hourly photospheric magnetograms using the Ni I 676.8 nm spectral line with an
    array of 242×256 pixels covering the solar disk. These magnetograms are used to derive
    synoptic maps which show a full-surface picture of the solar magnetic field.

    Notes
    -----
    SunPy automatically fixes non-compliant FITS metadata in GONG maps. Namely, CDELT2 is
    converted to degrees and CUNIT1 and CUNIT2 are set to 'deg'; DATE-OBS is converted to
    ISO format; BUNIT is set to 'G'; and HGLT_OBS, HGLN_OBS, and DSUN_OBS are set to the
    appropriate values for an observer on Earth.

    References
    ----------
    * `GONG Page <https://gong.nso.edu/>`_
    * `Magnetogram Synoptic Map Images Page <https://gong.nso.edu/data/magmap/>`_
    * `FITS header keywords <https://gong.nso.edu/data/DMAC_documentation/General/fitsdesc.html>`_
    * `Instrument Paper (pp. 203–208) <https://inis.iaea.org/collection/NCLCollectionStore/_Public/20/062/20062491.pdf>`_
    * `GONG+ Documentation <https://gong.nso.edu/data/DMAC_documentation/PipelineMap/GlobalMap.html>`_
    """

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        return (str(header.get('TELESCOP', '')).endswith('GONG') and
                str(header.get('CTYPE1', '').startswith('CRLN')))

    @property
    def date(self):
        return Time(self.meta.get('date-obs') + ' '
                    + self.meta.get('time-obs'))

    @property
    def scale(self):
        # Instead of the spacing in sin(lat), this should be 180/pi times
        # that value (see Thompson 2005)
        return SpatialPair(self.meta['cdelt1'] * self.spatial_units[0] / u.pixel,
                           self.meta['cdelt2'] * 180 / np.pi * self.spatial_units[0] / u.pixel)

    @property
    def unit(self):
        unit_str = self.meta.get('bunit', None)
        if unit_str is None:
            return
        elif unit_str == 'Gauss':
            return u.G
        return u.Unit(unit_str)

    @property
    def spatial_units(self):
        return SpatialPair(u.deg, u.deg)

    @property
    def observer_coordinate(self):
        observer_coord = get_earth(self.date)
        new_obs_frame = HeliographicStonyhurst(obstime=self.date)
        observer_coord = observer_coord.transform_to(new_obs_frame)

        sc = SkyCoord(
            obstime=self.date,
            lon=observer_coord.lon.to_value(u.deg),
            lat=observer_coord.lat.to_value(u.deg),
            radius=observer_coord.radius.to_value(u.m),
            unit=(u.deg, u.deg, u.m),
            frame="heliographic_stonyhurst"
        )
        sc = sc.heliographic_stonyhurst

        return SkyCoord(sc.replicate(rsun=self._rsun_meters(sc.radius)))
