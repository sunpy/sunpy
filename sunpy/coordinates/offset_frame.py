import astropy.units as u
from astropy.coordinates import SkyOffsetFrame, SphericalRepresentation

__all__ = ['NorthOffsetFrame']


class NorthOffsetFrame:
    """
    A frame which is offset from another frame such that it shares the same origin, but has its
    "north pole" (i.e., the Z axis) in a different direction.

    The original coordinate frame and the direction of the new north pole are specified by the
    ``north`` keyword.

    This class should be used when specifying a new north pole is natural.  In constrast, for
    shifting the origin in the projected sky (e.g., where helioprojective X and Y coordinates are
    zero), use `~astropy.coordinates.SkyOffsetFrame` instead.

    Parameters
    ----------
    north : `~sunpy.coordinates.frames.HeliographicStonyhurst`
        The direction and frame for the new "north pole".

    Examples
    --------
    .. minigallery:: sunpy.coordinates.NorthOffsetFrame

    Notes
    -----
    ``NorthOffsetFrame`` is a wrapper around the
    `~astropy.coordinates.SkyOffsetFrame` factory class. This class will
    calculate the desired coordinates of the ``origin`` from the ``north``
    keyword argument and then create a `~astropy.coordinates.SkyOffsetFrame`.

    Using this frame is equivalent to using
    `~astropy.coordinates.SkyOffsetFrame` with ``lat = lat - 90*u.deg`` for a
    position of the north pole in the original northern hemisphere.
    """

    def __new__(cls, *args, **kwargs):
        origin_frame = kwargs.pop('north', None)
        if origin_frame is None:
            raise TypeError("Can't initialize an NorthOffsetFrame without a `north` keyword.")
        if hasattr(origin_frame, 'frame'):
            origin_frame = origin_frame.frame

        rep = origin_frame.spherical

        lon = rep.lon
        lat = rep.lat
        if lat > 0*u.deg:
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
