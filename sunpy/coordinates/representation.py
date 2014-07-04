"""
SunPy's custom representation classes.
Part of the proposed Coordinates API.
@author: Pritish C. (VaticanCameos)
"""

# NumPy
import numpy as np

# Astropy
from astropy.extern import six
from astropy.utils.compat.odict import OrderedDict
from astropy import units as u
from astropy.coordinates.representation import SphericalRepresentation
from astropy.coordinates import Longitude, Latitude, Distance

class SphericalRepresentation180(SphericalRepresentation):
    """
    Representation of points in 3D Spherical coordinates.
    This representation allows for a negative Longitude.
    It does so by setting wrap_angle=180 degrees.

    Parameters
    ----------
    lon, lat: `~astropy.units.Quantity`
        The longitude and latitude of the point(s) in angular units. The
        latitude should be between -90 and +90 degrees, and the longitude
        is allowed to have any value between -360 to 360 degrees. These
        can also be instances of `~astropy.units.Angle`,
        `~astropy.coordinates.Longitude`, or `~astropy.coordinates.Latitude`.

    distance: `~astropy.units.Quantity`
        The distance to the point(s). If the distance is a length, it is
        passed to the :class:`~astropy.coordinates.Distance` class, else
        it is passed to the :class:`~astropy.units.Quantity` class.

    copy: bool, optional
        If True, arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('lon', Longitude),
                                ('lat', Latitude),
                                ('distance', u.Quantity)])
    recommended_units = {'lon': u.deg, 'lat': u.deg}

    def __init__(self, lon, lat, distance, copy=True):

        if not isinstance(lon, u.Quantity) or isinstance(lon, Latitude):
            raise TypeError('lon should be a Quantity, Angle or Longitude.')
        if not isinstance(lat, u.Quantity) or isinstance(lat, Longitude):
            raise TypeError('lat should be a Quantity, Angle or Latitude.')

        lon = Longitude(lon, copy=copy, wrap_angle=180*u.deg)
        lat = Latitude(lat, copy=copy)

        distance = u.Quantity(distance, copy=copy)
        if distance.unit.physical_type == 'length':
            distance = distance.view(Distance)

        try:
            lon, lat, distance = broadcast_quantity(lon, lat, distance, copy=copy)
        except:
            raise ValueError("Input parameters lon, lat and distance cannot be broadcast.")

        self._lon = lon
        self._lat = lat
        self._distance = distance

    @property
    def lon(self):
        """
        The longitude of the point(s).
        """
        return self._lon

    @property
    def lat(self):
        """
        The latitude of the point(s).
        """
        return self._lat

    @property
    def distance(self):
        """
        The distance to the point(s).
        """
        return self._distance
