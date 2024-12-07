"""
This file defines the net attributes that can be used to search the SOAR.
"""

import warnings

import astropy.units as u
import sunpy.net.attrs as a
from astropy.units import quantity_input
from sunpy.net.attr import (AttrAnd, AttrOr, AttrWalker, DataAttr, Range,
                            SimpleAttr)
from sunpy.util.exceptions import SunpyUserWarning

__all__ = ["SOOP", "Distance", "Product"]


class Product(SimpleAttr):
    """
    The data product descriptor to search for.

    Makes the value passed lower so that it is case insensitive as all
    descriptors on the SOAR are now lowercase.
    """

    def __init__(self, value) -> None:
        self.value = value.lower()


class SOOP(SimpleAttr):
    """
    The SOOP name to search for.
    """


class Distance(Range):
    """
    Specifies the distance range.

    Parameters
    ----------
    dist_min : `~astropy.units.Quantity`
        The lower bound of the range.
    dist_max : `~astropy.units.Quantity`
        The upper bound of the range.

    Notes
    -----
    The valid units for distance are AU, km, and mm. Any unit directly
    convertible to these units is valid input. This class filters the query
    by solar distance without relying on a specific distance column.
    """

    @quantity_input(dist_min=u.m, dist_max=u.m)
    def __init__(self, dist_min: u.Quantity, dist_max: u.Quantity):
        # Ensure both dist_min and dist_max are scalar values
        if not all([dist_min.isscalar, dist_max.isscalar]):
            msg = "Both dist_min and dist_max must be scalar values."
            raise ValueError(msg)

        target_unit = u.AU
        # Convert both dist_min and dist_max to the target unit
        dist_min = dist_min.to(target_unit)
        dist_max = dist_max.to(target_unit)

        super().__init__(dist_min, dist_max)

    def collides(self, other):
        """
        Check if the other attribute collides with this attribute.
        """
        return isinstance(other, self.__class__)


walker = AttrWalker()


@walker.add_creator(AttrOr)
def create_or(wlk, tree):
    """
    Creator for OR.

    Loops through the next level down in the tree and appends the
    individual results to a list.
    """
    return [wlk.create(sub) for sub in tree.attrs]


@walker.add_creator(AttrAnd, DataAttr)
def create_and(wlk, tree):
    """
    Creator for And and other simple attributes.

    No walking needs to be done, so simply call the applier function.
    """
    result = []
    wlk.apply(tree, result)
    return [result]


@walker.add_applier(AttrAnd)
def apply_and(wlk, and_attr, params) -> None:
    """
    Applier for And.

    Parameters
    ----------
    wlk : AttrWalker
    and_attr : AttrAnd
        The AND attribute being applied. The individual attributes being
        AND'ed together are accessible with ``and_attr.attrs``.
    params : list[str]
        List of search parameters.
    """
    for iattr in and_attr.attrs:
        wlk.apply(iattr, params)


"""
Below are appliers for individual attributes.

The all convert the attribute object into a query string, that will eventually
be passed as a query to the SOAR server. They all have the signature:

Parameters
----------
wlk : AttrWalker
    The attribute walker.
attr :
    The attribute being applied.
params : list[str]
    List of search parameters.
"""


@walker.add_applier(a.Time)
def _(wlk, attr, params) -> None:
    start = attr.start.strftime("%Y-%m-%d %H:%M:%S")
    end = attr.end.strftime("%Y-%m-%d %H:%M:%S")
    params.append(f"begin_time>='{start}' AND begin_time<='{end}'")


@walker.add_applier(a.Level)
def _(wlk, attr, params) -> None:
    level = attr.value
    if isinstance(level, int):
        level = f"L{level}"

    level = level.upper()
    allowed_levels = ("L0", "L1", "L2", "L3", "LL01", "LL02", "LL03")
    if level not in allowed_levels:
        warnings.warn(
            f"level not in list of allowed levels for SOAR: {allowed_levels}",
            SunpyUserWarning,
            stacklevel=2,
        )

    params.append(f"level='{level}'")


@walker.add_applier(a.Instrument)
def _(wlk, attr, params) -> None:
    params.append(f"instrument='{attr.value}'")


@walker.add_applier(Product)
def _(wlk, attr, params) -> None:
    params.append(f"descriptor='{attr.value}'")


@walker.add_applier(a.Provider)
def _(wlk, attr, params) -> None:
    params.append(f"provider='{attr.value}'")


@walker.add_applier(SOOP)
def _(wlk, attr, params) -> None:
    params.append(f"soop_name='{attr.value}'")


@walker.add_applier(a.Detector)
def _(wlk, attr, params) -> None:
    params.append(f"Detector='{attr.value}'")


@walker.add_applier(a.Wavelength)
def _(wlk, attr, params) -> None:
    wavemin = attr.min.value
    wavemax = attr.max.value
    params.append(f"Wavemin='{wavemin}' AND Wavemax='{wavemax}'")


@walker.add_applier(Distance)
def _(wlk, attr, params):
    # The `Distance` attribute is used to filter the query by solar distance
    # without relying on a specific distance column. It is commonly used
    # to filter the query without time consideration.
    dmin = attr.min.value
    dmax = attr.max.value
    min_possible = 0.28
    max_possible = 1.0

    if not (min_possible <= dmin <= max_possible) or not (min_possible <= dmax <= max_possible):
        warnings.warn(
            "Distance values must be within the range 0.28 AU to 1.0 AU.",
            SunpyUserWarning,
            stacklevel=2,
        )
    params.append(f"DISTANCE({dmin},{dmax})")
