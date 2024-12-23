"""
ADAPT Map subclass definitions
"""


from sunpy.map.mapbase import GenericMap, SpatialPair
from sunpy.time import parse_time

__all__ = ['ADAPTMap']

class ADAPTMap(GenericMap):
    """
    ADAPT magnetic flux map

    The ADAPT (Air Force Data Assimilative Photospheric Flux Transport) flux transport model evolves an ensemble of realizations, using flux transport processes during periods for which there are no observations, and updates the ensemble using the ensemble least-squares (EnLS) data assimilation method to account for model and observational uncertainties.
    These data products are Carrington maps with an additional third dimension corresponding to the ensemble of models used to produce the ADAPT map.

    Here, the third dimension is ignored and only the first element of the third dimension is used to create the  resulting map.

    .. warning::

        Please note that there are several data arrays within the ADAPT FITS files and `~sunpy.map.Map` will, by default, try to read them all and fail.
        In these cases, you must specify the header-data Pair you want to read.
        For these data, it will always be the first one.
        You can specify this by passing the ``hdus`` keyword argument to `~sunpy.map.Map`

        .. code-block:: python

            >>> sunpy.map.Map("adapt40311_03k012_202401020800_i00005600n1.fts.gz", hdus=0)  # doctest: +SKIP

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
        return self._get_date('date-obs') or self._get_date('maptime') or super().date

    def _set_date(self, date):
        self.meta['date-obs'] = parse_time(date).utc.isot
        if 'maptime' in self.meta:
            del self.meta['maptime']

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an ADAPT map."""
        return header.get('model') == 'ADAPT'
