"""
ADAPT Map subclass definitions
"""

from astropy.time import Time

from sunpy.map.mapbase import GenericMap, SpatialPair

__all__ = ['ADAPTMap']

class ADAPTMap(GenericMap):
    """
    ADAPT magnetic flux map

    The ADAPT (Air Force Data Assimilative Photospheric Flux Transport) flux transport model evolves an ensemble of realizations, using flux transport processes during periods for which there are no observations, and updates the ensemble using the ensemble least-squares (EnLS) data assimilation method to account for model and observational uncertainties.
    These data products are Carrington maps with an additional third dimension corresponding to the ensemble of models used to produce the ADAPT map.

    The third dimension is ignored and only the first element of the third dimension is used to create the map.

    References
    ----------
    * `NSO ADAPT Page <https://nso.edu/data/nisp-data/adapt-maps/>`__
    * `Data archive <https://gong.nso.edu/adapt/maps/>`__
    """
    @property
    def coordinate_system(self):
        ctype1, ctype2 = self.meta['ctype1'], self.meta['ctype2']
        if ctype1 == 'Long':
            ctype1 = 'CRLN-CAR'
        if ctype2 == 'Lat':
            ctype2 = 'CRLT-CAR'
        return SpatialPair(ctype1, ctype2)

    @property
    def date(self):
        return Time(self.meta.get('date-obs') or self.meta.get('maptime') or super().date)


    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an ADAPT map."""
        return header.get('model') == 'ADAPT'
