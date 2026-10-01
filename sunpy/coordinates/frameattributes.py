import datetime
from contextlib import contextmanager

import astropy.units as u
from astropy.coordinates import BaseCoordinateFrame, CoordinateAttribute, SkyCoord, TimeAttribute
from astropy.time import Time

from sunpy.time import parse_time

try:
    from astropy.coordinates import impose_frame_attributes
    _AssumedAttributeMixin = object()
except ImportError:
    from ._impose_attribute_backport import _AssumedAttributeMixin, impose_frame_attributes


__all__ = ["impose_observer", "TimeFrameAttributeSunPy", "ObserverCoordinateAttribute"]


@contextmanager
def impose_observer(observer):
    """
    Impose the provided observer and it's obstime on all frames with an observer attribute.

    Parameters
    ----------
    observer: astropy.coordinates.SkyCoord
        The observer to impose, should have an obstime which will also be imposed.

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> import sunpy.coordinates
    >>> from sunpy.coordinates import Helioprojective
    >>> import astropy.units as u
    >>> sc = SkyCoord(0*u.deg, 0*u.deg, 5*u.km,
    ...               obstime="2010/01/01T00:00:00", observer="earth", frame="helioprojective")
    >>> sc
    <SkyCoord (Helioprojective: obstime=2010-01-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, km)
        (0., 0., 5.)>

    This transformation does an observer shift and round trip through Heliographic coordinates:

    >>> target_frame = Helioprojective(observer="mars", obstime=sc.obstime)
    >>> sc.transform_to(target_frame)
    <SkyCoord (Helioprojective: obstime=2010-01-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'mars'>): (Tx, Ty, distance) in (arcsec, arcsec, km)
        (-79535.05476998, -424.5843386, 1.105226e+08)>

    This transformation does not, because the observer of the input and output coordinates are the same, so in this case the transform does nothing.

    >>> with impose_observer(target_frame.observer):
    ...     sc.transform_to(Helioprojective(obstime="2010/01/02T00:00:00"))
    <SkyCoord (Helioprojective: obstime=2010-01-01T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'mars'>): (Tx, Ty, distance) in (arcsec, arcsec, km)
        (0., 0., 5.)>
    """
    with impose_frame_attributes(observer=observer, obstime=observer.obstime):
        yield

class TimeFrameAttributeSunPy(_AssumedAttributeMixin, TimeAttribute):
    """
    Frame attribute descriptor for quantities that are Time objects.
    See the `~astropy.coordinates.Attribute` API doc for further
    information.

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value
        if ``default`` is ``None`` and no value was supplied during initialization.

    Returns
    -------
    frame_attr : descriptor
        A new data descriptor to hold a frame attribute
    """

    def convert_input(self, value):
        """
        Convert input value to a Time object and validate by running through the
        Time constructor. Also check that the input was a scalar.

        Parameters
        ----------
        value : object
            Input value to be converted.

        Returns
        -------
        out, converted : correctly-typed object, boolean
            Tuple consisting of the correctly-typed object and a boolean which
            indicates if conversion was actually performed.

        Raises
        ------
        ValueError
            If the input is not valid for this attribute.
        """
        if value is None:
            return None, False

        elif isinstance(value, Time):
            out = value
            converted = False

        elif isinstance(value, str):
            if value == 'now':
                return Time(datetime.datetime.now()), True

            try:
                out = Time(parse_time(value))
            except Exception as err:
                raise ValueError(f'Invalid time input {self.name}={value!r}\n{err}')
            converted = True
        else:
            try:
                out = Time(value)
            except Exception as err:
                raise ValueError(f'Invalid time input {self.name}={value!r}\n{err}')
            converted = True

        return out, converted


class ObserverCoordinateAttribute(_AssumedAttributeMixin, CoordinateAttribute):
    """
    An Attribute to describe the location of the observer in the solar system.
    The observer location can be given as a string of a known observer, which
    will be converted to a coordinate as long as the ``obstime`` attribute is
    valid on the instance of the frame. Alternatively a low-level frame class
    *or* a `~astropy.coordinates.SkyCoord` can be provided to specify the
    location of the observer. If a `~astropy.coordinates.SkyCoord` is passed it
    will always be converted to the low-level frame class when accessed.

    Parameters
    ----------
    frame : a coordinate frame class
        The type of frame this attribute can be
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    """

    def convert_input(self, value):
        # Keep string here.
        if isinstance(value, str):
            return value, False
        else:
            # Make sure that the coordinate is 3D
            if hasattr(value, 'make_3d'):
                value = value.make_3d()

            # Upgrade the coordinate to a `SkyCoord` so that frame attributes will be merged
            if isinstance(value, BaseCoordinateFrame) and not isinstance(value, self._frame):
                value = SkyCoord(value)

            result, converted = super().convert_input(value)

            # If already in the correct frame, make sure it has the default representation
            if (not converted and isinstance(result, BaseCoordinateFrame) and
                    not issubclass(result.representation_type, result.default_representation)):
                result = result.replicate(representation_type=result.default_representation)
                converted = True

            return result, converted

    def _convert_string_to_coord(self, out, obstime):
        """
        Given a value and and frame instance calculate the position of the
        object given as a string.
        """

        # Import here to prevent circular import
        from .ephemeris import get_body_heliographic_stonyhurst

        obscoord = get_body_heliographic_stonyhurst(out, obstime)

        if out == "earth":
            rep = obscoord.spherical
            rep.lon[()] = 0*u.deg
            obscoord = obscoord.realize_frame(rep)

        return obscoord

    def __get__(self, instance, frame_cls=None):
        # If instance is None then we can't get obstime so it doesn't matter.
        if instance is not None:
            observer = getattr(instance, '_' + self.name)
            obstime = getattr(instance, 'obstime', None)  # TODO: Why is this `None` needed?

            # If the observer is a string and we have obstime then calculate
            # the position of the observer.
            if isinstance(observer, str):
                if observer != "self" and obstime is not None:
                    new_observer = self._convert_string_to_coord(observer.lower(), obstime)
                    new_observer.object_name = observer
                    setattr(instance, '_' + self.name, new_observer)
                else:
                    return observer

        return super().__get__(instance, frame_cls=frame_cls)
