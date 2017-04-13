import astropy.units as u
from astropy.coordinates import SkyOffsetFrame, SphericalRepresentation, UnitSphericalRepresentation

__all__ = ['NorthOffsetFrame']


class NorthOffsetFrame(object):
    """
    A frame which is relative to some position and another frame. Based on
    `astropy.coordinates.SkyOffsetFrame`

    Coordinates in a NorthOffsetFrame are both centered on the position
    specified by the ``north`` keyword *and* they are oriented in the same
    manner as the ``north`` frame.

    Unlike `~astropy.coordinates.SkyOffsetFrame` a `NorthOffsetFrame` allows
    you to specify the position of the new north pole rather than the new
    origin to centre the new frame.

    Examples
    --------

    A common use for this is to create a frame derived from Heliographic, which
    has the north pole at some point of interest. In this new frame, lines of
    longitude form great circles radially away from the point, and lines of
    latitude measure angular distance from the point.

    In this example the new frame is shifted so the new north pole is at (20,
    20) in the Heliographic Stonyhurst frame. The new grid is overplotted in
    blue.

    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        import astropy.units as u

        from astropy.coordinates import SkyCoord
        from sunpy.coordinates import NorthOffsetFrame

        import sunpy.map
        from sunpy.data.sample import AIA_171_IMAGE

        m = sunpy.map.Map(AIA_171_IMAGE)

        north = SkyCoord(20*u.deg, 20*u.deg, frame="heliographic_stonyhurst")
        new_frame = NorthOffsetFrame(north=north)

        ax = plt.subplot(projection=m)
        m.plot()

        overlay = ax.get_coords_overlay('heliographic_stonyhurst')
        overlay[0].set_ticks(spacing=30. * u.deg, color='white')
        overlay.grid(ls='-', color='white')

        overlay = ax.get_coords_overlay(new_frame)
        overlay[0].set_ticks(spacing=30. * u.deg)
        overlay.grid(ls='-', color='blue')

    Notes
    -----
    ``NorthOffsetFrame`` is a wrapper around the
    `~astropy.coordinates.SkyOffsetFrame` factory class. This class will
    calculate the desired coordianates of the ``origin`` from the ``north``
    keyword argument and then create a `~astropy.coordinates.SkyOffsetFrame`.

    Using this frame is equivalent to using
    `~astropy.coordinates.SkyOffsetFrame` with ``lat = lat - 90*u.deg`` for a
    position of the north pole in the original northern hemisphere.

    This class will only work for Heliographic-Stonyhurst and Heliographic
    Carrington frames, and not helioprojective. If you want to rotate a
    helioprojective frame, it would be more natural to use the
    `~astropy.coordinates.SkyOffsetFrame`.
    """

    def __new__(cls, *args, **kwargs):
        origin_frame = kwargs.pop('north', None)
        if origin_frame is None:
            raise TypeError("Can't initialize an NorthOffsetFrame without origin= keyword.")
        if hasattr(origin_frame, 'frame'):
            origin_frame = origin_frame.frame

        if not isinstance(origin_frame.data, SphericalRepresentation):
            rep = origin_frame.represent_as(SphericalRepresentation)
        else:
            rep = origin_frame.data

        lon = rep.lon
        lat = rep.lat
        if lat > 0*u.deg:
            lat = lat - 90*u.deg
            rotation = None
        else:
            lon = lon - 180*u.deg
            lat = -90*u.deg - lat
            rotation = 180*u.deg

        if isinstance(origin_frame.data, UnitSphericalRepresentation):
            new_rep = origin_frame.representation(lon=lon,
                                                  lat=lat)

        else:
            new_rep = origin_frame.representation(lon=lon,
                                                  lat=lat,
                                                  distance=rep.distance)

        new_origin = origin_frame.realize_frame(new_rep)
        kwargs['origin'] = new_origin
        kwargs['rotation'] = rotation

        return SkyOffsetFrame(*args, **kwargs)
