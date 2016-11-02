import astropy.units as u
from astropy.coordinates import SkyOffsetFrame, SphericalRepresentation

__all__ = ['NorthOffsetFrame']


class NorthOffsetFrame():
    """
    A frame which is relative to some position and another frame. Based on
    `astropy.coordinates.SkyOffsetFrame`

    Coordinates in a NorthOffsetFrame are both centered on the position
    specified bu the ``north`` keyword *and* they are oriented in the same
    manner as the ``north`` frame.

    Unlike `~astropy.coordinates.SkyOffsetFrame` a `NorthOffsetFrame` allows
    you to specify the position of the new north pole rather than the new
    origin to centre the new frame.

    Notes
    -----
    ``NorthOffsetFrame`` is a wrapper around the
    `~astropy.coordinates.SkyOffsetFrame` factory class. This class will
    calculate the desired coordianates of the ``origin`` from the ``north``
    keyword argument and then create a `~astropy.coordinates.SkyOffsetFrame`.

    Using this frame is equivalent to using
    `~astropy.coordinates.SkyOffsetFrame` with ``lat = lat - 90*u.deg``
    """

    def __new__(cls, *args, **kwargs):
        origin_frame = kwargs.pop('north', None)
        if origin_frame is None:
            raise TypeError("Can't initialize an NorthOffsetFrame without origin= keyword.")
        if hasattr(origin_frame, 'frame'):
            origin_frame = origin_frame.frame

        rep = origin_frame.represent_as(SphericalRepresentation)
        lon = rep.lon
        lat = rep.lat
        if rep.lat > 0*u.deg:
            lat = lat - 90*u.deg
            rotation = None
        else:
            lon = lon - 180*u.deg
            lat = -90*u.deg - lat
            rotation = 180*u.deg
        new_rep = SphericalRepresentation(lon=lon,
                                          lat=lat,
                                          distance=rep.distance)

        new_origin = origin_frame.realize_frame(new_rep)
        kwargs['origin'] = new_origin
        kwargs['rotation'] = rotation
        return SkyOffsetFrame(*args, **kwargs)
