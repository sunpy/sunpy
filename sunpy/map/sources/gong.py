"""
GONG Map subclass definitions
"""
import numpy as np

import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord

from sunpy.coordinates import get_earth
from sunpy.map import GenericMap
from sunpy.time import parse_time

__all__ = ['GONGSynopticMap', 'GONGHalphaMap']

from sunpy.map.mapbase import SpatialPair

_SITE_NAMES = {
    'LE': 'Learmonth',
    'UD': 'Udaipur',
    'TD': 'El Teide',
    'CT': 'Cerro Tololo',
    'BB': 'Big Bear',
    'ML': 'Mauna Loa'
}


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
    If you have ``pfsspy`` installed this map source will be used instead of the one built into ``pfsspy``.
    References
    ----------
    * `GONG Page <https://gong.nso.edu/>`__
    * `Magnetogram Synoptic Map Images Page <https://gong.nso.edu/data/magmap/>`__
    * `FITS header keywords <https://gong.nso.edu/data/DMAC_documentation/General/fitsdesc.html>`__
    * `Instrument Paper (pp. 203–208) <https://inis.iaea.org/collection/NCLCollectionStore/_Public/20/062/20062491.pdf>`__
    * `GONG+ Documentation <https://gong.nso.edu/data/DMAC_documentation/PipelineMap/GlobalMap.html>`__
    """

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        return (str(header.get('TELESCOP', '')).endswith('GONG') and
                str(header.get('CTYPE1', '')).startswith('CRLN'))

    @property
    def date(self):
        # The FITS file has a date that is made from the date-obs and time-obs keywords
        # Which is not what date-obs is supposed to be.
        if date_obs := self.meta.get('date-obs'):
            if time_obs := self.meta.get('time-obs'):
                date_obs = f"{date_obs} {time_obs}"
            date_obs = parse_time(date_obs)
        return date_obs or super().date

    def _set_date(self, date):
        if 'time-obs' in self.meta:
            self.meta['date-obs'], self.meta['time-obs'] = parse_time(date).utc.isot.split('T')
        else:
            self.meta['date-obs'] = parse_time(date).utc.isot

    @property
    def reference_date(self):
        return self.date

    def _set_reference_date(self, date):
        self._set_date(date)

    @property
    def scale(self):
        # Since, this map uses the cylindrical equal-area (CEA) projection,
        # the spacing should be modified to 180/pi times the original value
        # Reference: Section 5.5, Thompson 2006
        return SpatialPair(self.meta['cdelt1'] * self.spatial_units[0] / u.pixel,
                           self.meta['cdelt2'] * 180 / np.pi * self.spatial_units[0] / u.pixel)

    @property
    def spatial_units(self):
        return SpatialPair(u.deg, u.deg)

    @property
    def observer_coordinate(self):
        return get_earth(self.date)


class GONGHalphaMap(GenericMap):
    """
    GONG H-Alpha Map.

    The Global Oscillation Network Group (GONG) operates a six-station network of H-alpha
    imagers located around the Earth that observe the Sun nearly continuously.

    References
    ----------
    * `GONG H-Alpha Page <https://nso.edu/data/nisp-data/h-alpha/>`__
    * `GONG H-Alpha Observation Details <https://nispdata.nso.edu/webProdDesc2/presenter.php?file=halpha_fulldisk_images_overview.html&echoExact=0&name=Overview%20:%20GONG%20H-alpha%20Full-disk%20Images>`__
    * `GONG Header Keywords <https://gong.nso.edu/data/HEADER_KEY.html>`__
    * `GONG H-alpha Full-disk Images. <https://doi.org/10.25668/as28-7p13>`__
    """

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        return (str(header.get('TELESCOP', '')).endswith('GONG') and
                str(header.get('IMTYPE', '')).startswith('H-ALPHA'))

    @property
    def scale(self):
        solar_r = self.meta['SOLAR-R'] * u.arcsec
        return SpatialPair(solar_r / (self.meta['FNDLMBMI'] * u.pixel),
                           solar_r/ (self.meta['FNDLMBMA'] * u.pixel))

    @property
    def rsun_obs(self):
        # Header contains a radius keyword which seems to have a higher priority but for GONG Ha is in pixels
        return self.meta['SOLAR-R'] * u.arcsec

    @property
    def coordinate_system(self):
        """
        Coordinate system used

        Overrides the values in the header which are not understood by Astropy WCS
        """
        return SpatialPair("HPLN-TAN", "HPLT-TAN")

    @property
    def nickname(self):
        site = _SITE_NAMES.get(self.meta.get("sitename", ""), "UNKNOWN")
        return f'{self.observatory}, {site}'

    @property
    def spatial_units(self):
        return SpatialPair(u.deg, u.deg)

    @property
    def _earth_location(self):
        """Location of the observatory on Earth"""
        return EarthLocation.from_geodetic(lat=self.meta['site-lat'] * u.deg, lon=self.meta['site-lon'] * u.deg)

    @property
    def observer_coordinate(self):
        return SkyCoord(self._earth_location.get_itrs(self.date)).heliographic_stonyhurst
