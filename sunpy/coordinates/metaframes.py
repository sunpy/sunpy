"""
Coordinate frames that are defined relative to other frames
"""

import astropy.units as u
from astropy.coordinates.attributes import Attribute, QuantityAttribute
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import FunctionTransform

from sunpy.time import parse_time
from sunpy.time.time import _variables_for_parse_time_docstring
from sunpy.util.decorators import add_common_docstring
from ._transformations import _transformation_debug
from .frames import HeliographicStonyhurst, SunPyBaseCoordinateFrame
from .offset_frame import NorthOffsetFrame

__all__ = ['NorthOffsetFrame', 'RotatedSunFrame']


# The code for NorthOffsetFrame currently lives in `offset_frame.py`
# This changes its module to this file so that the docs think that the class is local
NorthOffsetFrame.__module__ = __name__


_rotatedsun_cache = {}


def _make_rotatedsun_cls(framecls):
    """
    Create a new class that is the rotated-Sun frame for a specific class of
    base frame. If such a class has already been created for this frame, the
    same class will be returned.

    This function is necessary because frame transformations depend
    on connection between specific frame *classes*. So each type of frame
    needs its own distinct rotated-Sun frame class. This function generates
    just that class, as well as ensuring that only one example of such a class
    actually gets created in any given Python session.
    """
    # This code reuses significant code from Astropy's implementation of SkyOffsetFrame
    # See licenses/ASTROPY.rst

    if framecls in _rotatedsun_cache:
        return _rotatedsun_cache[framecls]

    members = {'__doc__': f'{RotatedSunFrame.__doc__}\n{framecls.__doc__}'}

    # Copy over the defaults from the input frame class
    attrs_to_copy = ['_default_representation',
                     '_default_differential',
                     '_frame_specific_representation_info']
    for attr in attrs_to_copy:
        members[attr] = getattr(framecls, attr)

    _RotatedSunFramecls = type(f'RotatedSun{framecls.__name__}', (RotatedSunFrame,), members)

    @frame_transform_graph.transform(FunctionTransform, _RotatedSunFramecls, _RotatedSunFramecls)
    @_transformation_debug(f"{_RotatedSunFramecls.__name__}->{_RotatedSunFramecls.__name__}")
    def rotatedsun_to_rotatedsun(from_rotatedsun_coord, to_rotatedsun_frame):
        """Transform between two rotated-Sun frames."""
        # This transform goes through the parent frames on each side.
        # from_frame -> HGS -> to_frame
        int_frame = HeliographicStonyhurst(obstime=from_rotatedsun_coord.base.obstime)
        int_coord = from_rotatedsun_coord.transform_to(int_frame)
        return int_coord.transform_to(to_rotatedsun_frame)

    @frame_transform_graph.transform(FunctionTransform, HeliographicStonyhurst, _RotatedSunFramecls)
    @_transformation_debug(f"HGS->{_RotatedSunFramecls.__name__}")
    def reference_to_rotatedsun(hgs_coord, rotatedsun_frame):
        int_frame = HeliographicStonyhurst(obstime=rotatedsun_frame.base.obstime)
        int_coord = hgs_coord.make_3d().transform_to(int_frame)  # obstime change handled here

        # Rotate the coordinate in HGS
        int_coord = int_coord._apply_diffrot(-rotatedsun_frame.duration,
                                             rotatedsun_frame.rotation_model)

        # Transform from HGS
        new_coord = int_coord.transform_to(rotatedsun_frame.base)
        return rotatedsun_frame.realize_frame(new_coord.data)

    @frame_transform_graph.transform(FunctionTransform, _RotatedSunFramecls, HeliographicStonyhurst)
    @_transformation_debug(f"{_RotatedSunFramecls.__name__}->HGS")
    def rotatedsun_to_reference(rotatedsun_coord, hgs_frame):
        # Transform to HGS
        from_coord = rotatedsun_coord.base.realize_frame(rotatedsun_coord.data)
        if hasattr(from_coord, 'make_3d'):
            from_coord = from_coord.make_3d()
        int_frame = HeliographicStonyhurst(obstime=rotatedsun_coord.base.obstime)
        int_coord = from_coord.transform_to(int_frame)

        # Rotate the coordinate in HGS
        int_coord = int_coord._apply_diffrot(rotatedsun_coord.duration,
                                             rotatedsun_coord.rotation_model)

        # Transform from HGS
        return int_coord.transform_to(hgs_frame)  # obstime change handled here

    _rotatedsun_cache[framecls] = _RotatedSunFramecls
    return _RotatedSunFramecls


@add_common_docstring(**_variables_for_parse_time_docstring())
class RotatedSunFrame(SunPyBaseCoordinateFrame):
    """
    A frame that applies solar rotation to a base coordinate frame.

    .. note::

        See :ref:`sunpy-topic-guide-coordinates-rotatedsunframe` for how to use this class.

    In essence, the coordinate axes of the frame are distorted by differential solar rotation.
    This allows using a coordinate representation at one time (at the ``obstime`` of the base
    coordinate frame) to point to a location at a different time that has been differentially
    rotated by the time difference (``duration``).

    Parameters
    ----------
    representation : `~astropy.coordinates.BaseRepresentation` or ``None``
        A representation object or ``None`` to have no data. Alternatively, use coordinate
        component keyword arguments, which depend on the base frame.
    base : `~astropy.coordinates.SkyCoord` or low-level coordinate object.
        The coordinate which specifies the base coordinate frame. The frame must be a SunPy frame.
    duration : `~astropy.units.Quantity` or `~astropy.time.TimeDelta`
        The duration of solar rotation (defaults to zero days).
    rotated_time : {parse_time_types}
        The time to rotate the Sun to. If provided, ``duration`` will be set to the difference
        between this time and the observation time in ``base``.
    rotation_model : `str`
        Accepted model names are ``'howard'`` (default), ``'snodgrass'``, ``'allen'``, and ``'rigid'``.
        See the documentation for :func:`~sunpy.physics.differential_rotation.diff_rot` for differences
        between these models.

    Notes
    -----
    ``RotatedSunFrame`` is a factory class. That is, the objects that it
    yields are *not* actually objects of class ``RotatedSunFrame``. Instead,
    distinct classes are created on-the-fly for whatever the frame class is
    of ``base``.
    """
    # This code reuses significant code from Astropy's implementation of SkyOffsetFrame
    # See licenses/ASTROPY.rst

    # We don't want to inherit the frame attributes of SunPyBaseCoordinateFrame (namely `obstime`)
    # Note that this does not work for Astropy 4.3+, so we need to manually remove it below
    _inherit_descriptors_ = False

    # Even though the frame attribute `base` is a coordinate frame, we use `Attribute` instead of
    # `CoordinateAttribute` because we are preserving the supplied frame rather than converting to
    # a common frame.
    base = Attribute()

    duration = QuantityAttribute(default=0*u.day)
    rotation_model = Attribute(default='howard')

    def __new__(cls, *args, **kwargs):
        # We don't want to call this method if we've already set up
        # a rotated-Sun frame for this class.
        if not (issubclass(cls, RotatedSunFrame) and cls is not RotatedSunFrame):
            # We get the base argument, and handle it here.
            base_frame = kwargs.get('base', None)
            if base_frame is None:
                raise TypeError("Can't initialize a RotatedSunFrame without a `base` keyword.")

            # If a SkyCoord is provided, use the underlying frame
            if hasattr(base_frame, 'frame'):
                base_frame = base_frame.frame

            newcls = _make_rotatedsun_cls(base_frame.__class__)
            return newcls.__new__(newcls, *args, **kwargs)

        # http://stackoverflow.com/questions/19277399/why-does-object-new-work-differently-in-these-three-cases
        # See above for why this is necessary. Basically, because some child
        # may override __new__, we must override it here to never pass
        # arguments to the object.__new__ method.
        if super().__new__ is object.__new__:
            return super().__new__(cls)
        return super().__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        # Validate inputs
        base_frame = kwargs.get('base')
        if base_frame is None or base_frame.obstime is None:
            raise ValueError("The base coordinate frame must have a defined `obstime`.")

        duration = kwargs.pop('duration', None)
        rotated_time = kwargs.pop('rotated_time', None)

        if duration is not None and rotated_time is not None:
            raise ValueError("Specify either `duration` or `rotated_time`, not both.")

        if duration is not None:
            rotated_time = base_frame.obstime + duration

        if rotated_time is not None:
            kwargs['duration'] = (parse_time(rotated_time) - base_frame.obstime).to('day')

        super().__init__(*args, **kwargs)

        # Move data out from the base frame
        if self.base.has_data:
            if not self.has_data:
                # If the duration is an array but the data is scalar, upgrade data to an array
                if self.base.data.isscalar and not self.duration.isscalar:
                    self._data = self.base.data._apply('repeat', self.duration.shape)
                else:
                    self._data = self.base.data
            self._base = self.base.replicate_without_data()

    def as_base(self):
        """
        Returns a coordinate with the current representation and in the base coordinate frame.

        This method can be thought of as "removing" the
        `~sunpy.coordinates.metaframes.RotatedSunFrame` layer. Be aware that this method is not
        merely a coordinate transformation, because this method changes the location in inertial
        space that is being pointed to.
        """
        return self.base.realize_frame(self.data)

    @property
    def rotated_time(self):
        """
        Returns the sum of the base frame's observation time and the rotation of duration.
        """
        return self.base.obstime + self.duration

    def __reduce__(self):
        return (_rotatedsunframe_reducer, (self.base,), self.__dict__)


def _rotatedsunframe_reducer(base):
    return RotatedSunFrame.__new__(RotatedSunFrame, base=base)


# For Astropy 4.3+, we need to manually remove the `obstime` frame attribute from RotatedSunFrame
if 'obstime' in RotatedSunFrame.frame_attributes:
    del RotatedSunFrame.frame_attributes['obstime']
