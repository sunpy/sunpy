"""
GONG Map subclass definitions
"""
import numpy as np

import astropy.units as u
from astropy.time import Time

from sunpy.coordinates import get_earth
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
    If you have ``pfsspy`` installed this map source will be used instead of the one built into ``pfsspy``.

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
        return Time(f"{self.meta.get('date-obs')} {self.meta.get('time-obs')}")

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
