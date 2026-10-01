# TODO: Drop this whole file when min astropy is 8.1.0
from collections import defaultdict
from contextlib import contextmanager
from contextvars import ContextVar

import numpy as np

from astropy.coordinates import Attribute
from astropy.utils import ShapedLikeNDArray

# Skip doctests as it uses dev astropy
__doctest_skip__ = ["impose_frame_attributes"]


_imposed_attributes: ContextVar[defaultdict | None] = ContextVar(
    "_imposed_attributes", default=None
)


def _get_imposed_attributes():
    """A helper function so that we don't have a mutable default in ContextVar."""
    return _imposed_attributes.get() or {}



@contextmanager
def impose_frame_attributes(**kwargs):
    """
    Assume a value of one or more frame attributes.

    Parameters
    ----------
    kwargs
        This function accepts any keyword arguments and will override
        any frame attribute with a name matching that of the keyword
        name.

    Examples
    --------
    .. doctest-skip-all::


    >>> import astropy.units as u
    >>> from astropy.coordinates import EarthLocation, SkyCoord
    >>> from astropy.coordinates.attributes import impose_frame_attributes

    >>> location = EarthLocation(-5466045*u.m, -2404388*u.m, 2242133*u.m)

    Say you have two coordinates that represent points in images which are very close in time but not exactly the same:

    >>> coord1 = SkyCoord(10, 20, unit=u.deg, obstime="2026-01-01T00:00:00", location=location, frame="altaz")
    >>> coord2 = SkyCoord(10, 20, unit=u.deg, obstime="2026-01-01T00:01:00", location=location, frame="altaz")

    If you transform these to other frames, they will give slightly different results:

    >>> coord1.transform_to("icrs")
    <SkyCoord (ICRS): (ra, dec) in deg
        (36.10412254, 80.47633131)>
    >>> coord2.transform_to("icrs")
    <SkyCoord (ICRS): (ra, dec) in deg
        (36.35158456, 80.47671193)>

    If however, you are happy to ignore this small change in observer
    location, and you want to avoid more complex transformations or
    other effects through the transform graph, you can do this:

    >>> with impose_frame_attributes(obstime="2026-01-01T00:00:00"):
    ...     coord2.transform_to("icrs")
    <SkyCoord (ICRS): (ra, dec) in deg
        (36.10412254, 80.47633131)>

    Another example is computing the separation between two `~.HCRS` frames with slightly different obstimes.

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord

    >>> coord1 = SkyCoord(10*u.deg, 20*u.deg, frame='hcrs', obstime='2026-01-01 00:00:00')
    >>> coord2 = SkyCoord(20*u.deg, 30*u.deg, frame='hcrs', obstime='2026-01-01 00:00:00.001')

    When you try and compute the separation an error is raised because of the origin shift.

    >>> print(coord1.separation(coord2))  # doctest: +SKIP
    ...
    astropy.units.errors.UnitsError: The input HCRS coordinates do not have length units. This probably means you created coordinates with lat/lon but no distance.  Heliocentric<->ICRS transforms cannot function in this case because there is an origin shift.

    However, using this context manager:

    >>> with impose_frame_attributes(obstime=coord1.obstime):
    ...     print(coord1.separation(coord2))
    13d28m54.20360928s

    Both of these examples are trivial when constructing coordinates
    like this, but when operating on a large number of coordinates
    read from a file or with coordinate transforms in packages such as
    ``reproject``, it becomes harder to make these changes manually.
    """
    # Get the current imposed attributes
    current_attrs = _get_imposed_attributes()
    # Build the new set, overriding current with kwargs
    new_attrs = {**current_attrs, **kwargs}
    # Set the new imposed attrs, but make sure it's still a defaultdict
    token = _imposed_attributes.set(defaultdict(lambda: None, new_attrs))
    try:
        yield
    finally:
        # Reset to previous value using token
        _imposed_attributes.reset(token)


class _AssumedAttributeMixin(Attribute):  # Inherit for type checking mainly
    # This class re-implements the astropy.coordinates.attributes.Attribute.__get__ method
    # to support assumed attributes, maybe this will be upstreamed at somepoint
    def __get__(self, instance, frame_cls=None):
        if instance is None:
            # Return the descriptor instance to enable the retrieval of the docstring
            return self

        if (imposed_value := _get_imposed_attributes().get(self.name)) is not None:
            out = imposed_value
        else:
            out = getattr(instance, "_" + self.name, self.default)

        if out is None:
            out = getattr(instance, self.secondary_attribute, self.default)

        out, converted = self.convert_input(out)
        if instance is not None:
            # None if instance (frame) has no data!
            instance_shape = getattr(instance, "shape", None)
            if instance_shape is not None and (
                getattr(out, "shape", ()) and out.shape != instance_shape
            ):
                # If the shapes do not match, try broadcasting.
                try:
                    if isinstance(out, ShapedLikeNDArray):
                        out = out._apply(
                            np.broadcast_to, shape=instance_shape, subok=True
                        )
                    else:
                        out = np.broadcast_to(out, instance_shape, subok=True)
                except ValueError:
                    # raise more informative exception.
                    raise ValueError(
                        f"attribute {self.name} should be scalar or have shape"
                        f" {instance_shape}, but it has shape {out.shape} and could not"
                        " be broadcast."
                    )

                converted = True

            if converted and not imposed_value:
                setattr(instance, "_" + self.name, out)

        return out
