import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.coordinates.frames import Helioprojective
from sunpy.coordinates.screens import SphericalScreen
from sunpy.coordinates.utils import get_rectangle_coordinates
from sunpy.net._attrs import Time, Wavelength
from sunpy.net.attr import AttrAnd, AttrComparison, AttrOr, AttrWalker, DataAttr, SimpleAttr
from sunpy.util.exceptions import warn_user

__all__ = ['Series', 'Protocol', 'Notify', 'Segment', 'PrimeKey', 'Cutout', "Keyword"]


# Define a custom __dir__ to restrict tab-completion to __all__
def __dir__():
    return __all__


class Series(SimpleAttr):
    """
    The JSOC Series to Download.

    This is the list of `Series <http://jsoc.stanford.edu/JsocSeries_DataProducts_map.html>`__.
    """


class PrimeKey(DataAttr):
    """
    Prime Keys

    Parameters
    ----------
    label : str
    value : str
    """

    def __init__(self, label, value):
        super().__init__()
        self.label = label
        self.value = value

    def __repr__(self):
        return f"{object.__repr__(self)}" + "\n" + f"{self.label, self.value}"

    def collides(self, other):
        return False


class KeywordComparison(AttrComparison):
    """
    Allows comparison filtering of the JSOC Keywords with the ability to specify the comparison operator.

    Parameters
    ----------
    name : str
    operator : str
    value : Numeric
    """


class Keyword(SimpleAttr):
    """
    Allows comparison filtering of the JSOC Keywords.

    Parameters
    ----------
    value : str
    """

    def __lt__(self, other):
        return KeywordComparison(self.value, '<', other)

    def __le__(self, other):
        return KeywordComparison(self.value, '<=', other)

    def __gt__(self, other):
        return KeywordComparison(self.value, '>', other)

    def __ge__(self, other):
        return KeywordComparison(self.value, '>=', other)

    def __eq__(self, other):
        return KeywordComparison(self.value, '=', other)

    def __ne__(self, other):
        return KeywordComparison(self.value, '!=', other)

    def collides(self, other):
        return isinstance(other, Keyword)


class Segment(SimpleAttr):
    """
    Segments choose which files to download when there are more than
    one present for each record e.g. 'image'.
    """

    def collides(self, other):
        return False


class Protocol(SimpleAttr):
    """
    The type of download to request one of
    ("FITS", "JPEG", "MPG", "MP4", or "as-is").
    Only FITS is supported, the others will require extra keywords.
    """


class Notify(SimpleAttr):
    """
    An email address to get a notification to when JSOC has staged your request.
    """

    def __init__(self, value):
        super().__init__(value)
        if value is None:
            raise ValueError("Notify attribute must contain an email address")
        if value.find('@') == -1:
            raise ValueError("Notify attribute must contain an '@' symbol "
                             "to be a valid email address")
        self.value = value


class Cutout(DataAttr):
    """
    Select a cutout region.

    The JSOC allows for users to request cutouts. This process is performed server
    side so as to allow users to download only the portions of the full-disk images
    they are interested in. For a detailed explanation of the routine
    used to perform these cutouts on the JSOC server, see
    http://jsoc.stanford.edu/doxygen_html/group__im__patch.html.

    Parameters
    ----------
    bottom_left : `~astropy.coordinates.SkyCoord`
        Helioprojective coordinate for the bottom left corner of the cutout.
    top_right : `~astropy.coordinates.SkyCoord`, optional
        Coordinate for the top right corner of the cutout. If this is
        not specified, both ``width`` and ``height`` must both be specified.
    width : `~astropy.units.Quantity`, optional
        Width of the cutout. If this parameter, along with ``height``, is
        not specified, ``top_right`` must be specified.
    height : `~astropy.units.Quantity`, optional
        Height of the cutout. If this parameter, along with ``width``, is
        not specified, ``top_right`` must be specified.
    tracking : `bool`, optional
        If True, the field of view follows the rotation of the Sun
    register : `bool`, optional
        If True, use sub-pixel registration when cropping to the target location.
    nan_off_limb : `bool`, optional
        If True, all off-limb pixels are set to NaN

    See Also
    --------
    sunpy.coordinates.utils.get_rectangle_coordinates

    Notes
    -----
    The ``bottom_left`` coordinate must be in the `~sunpy.coordinates.Helioprojective`
    frame, and must have both ``obstime`` and ``observer`` defined. The supplied
    coordinates are transformed into the frame that the JSOC uses to interpret an
    ``arcsec`` cutout specification, rather than being assumed to already be in it.

    .. warning::

        The JSOC interprets the cutout specification as seen from SDO, but the
        transformation performed here uses **Earth** as the observer, as an
        approximation of SDO's location. SDO is in an inclined geosynchronous
        orbit, so it is displaced from Earth's center by at most one
        geosynchronous radius (~42,164 km). Because Helioprojective coordinates
        are referenced to the Sun-center direction, the bulk parallax from this
        displacement cancels, and only a depth-dependent residual remains: at
        most ~0.3 arcsec, or roughly half an AIA or HMI pixel. The residual is
        largest near disk center and falls to zero at the limb. If you need
        better accuracy than this, request a larger cutout and crop it locally.

    Transforming a two-dimensional Helioprojective coordinate between observers
    is ill-posed, because a 2D coordinate specifies only a line of sight rather
    than a point in space. An assumption about the distance to the coordinate is
    therefore required:

    * On-disk coordinates are assumed to lie on the solar surface, as defined by
      the ``rsun`` frame attribute. This is the usual `sunpy` assumption and is a
      good one for features in the photosphere/low corona.
    * Off-disk coordinates are placed on a `~sunpy.coordinates.SphericalScreen`
      centered on the original observer, i.e., they are assumed to lie at the
      same distance from the observer as the Sun's center.

    Because the corners of the requested rectangle do not in general remain a
    rectangle under this transformation, the returned cutout is the bounding box
    of the four transformed corners, so the requested field of view is fully
    contained in the result.

    A `~sunpy.util.exceptions.SunpyUserWarning` is emitted whenever a
    transformation is actually applied, since the result depends on the
    assumptions above. Coordinates that already have Earth as the observer are
    passed through unchanged and do not warn.

    If ``tracking`` is `True`, the center of the cutout is required to be on the
    solar disk, otherwise the JSOC will produce unexpected output.
    """
    @u.quantity_input
    def __init__(self, bottom_left, top_right=None, width: u.arcsec = None,
                 height: u.arcsec = None, tracking=False, register=False,
                 nan_off_limb=False):
        super().__init__()
        bl, tr = get_rectangle_coordinates(bottom_left, top_right=top_right, width=width, height=height)
        if not isinstance(bl.frame, Helioprojective):
            raise ValueError("`bottom_left` must be in the `Helioprojective` frame, but is instead "
                             f"in the `{bl.frame.__class__.__name__}` frame")
        if bl.obstime is None:
            raise ValueError("`bottom_left` must have `obstime` set, because the JSOC cutout "
                             "request requires a reference time.")
        if bl.frame.observer is None:
            raise ValueError("`bottom_left` must have `observer` set, because the coordinates "
                             "need to be transformed to a Helioprojective frame with Earth as "
                             "the observer. Use ``observer='earth'`` to reproduce the behavior "
                             "of sunpy <8.1, which assumed the coordinates needed no "
                             "transformation.")

        bl, tr = self._transform_to_earth_observer(bl, tr)

        center_x = (bl.Tx + tr.Tx) / 2
        center_y = (bl.Ty + tr.Ty) / 2
        center = SkyCoord(center_x, center_y, frame=bl.frame)
        if tracking:
            # import here so net won't depend on map
            from sunpy.map.maputils import coordinate_is_on_solar_disk
            if not coordinate_is_on_solar_disk(center):
                raise ValueError("Tracking is enabled, but the center of the cutout "
                                 f"(Tx={center_x}, Ty={center_y}) is not on the solar disk.")

        self.value = {
            't_ref': bl.obstime.isot,
            # JSOC input is disable tracking so take the negative
            't': int(not tracking),
            'r': int(register),
            'c': int(nan_off_limb),
            'locunits': 'arcsec',
            'boxunits': 'arcsec',
            'x': center_x.to_value('arcsec'),
            'y': center_y.to_value('arcsec'),
            'width': (tr.Tx - bl.Tx).to_value('arcsec'),
            'height': (tr.Ty - bl.Ty).to_value('arcsec'),
        }

    @staticmethod
    def _transform_to_earth_observer(bl, tr):
        """
        Transform the corners of the requested rectangle to a Helioprojective
        frame with Earth as the observer, which is used here as an approximation
        of SDO's location (see the class docstring for the size of the error).

        Returns the bottom-left and top-right corners of the bounding box of the
        four transformed corners. If the input observer is already Earth, the
        inputs are returned unchanged.
        """
        earth_frame = Helioprojective(obstime=bl.obstime, observer='earth', rsun=bl.frame.rsun)
        if bl.frame.is_equivalent_frame(earth_frame):
            return bl, tr

        corners = SkyCoord(u.Quantity([bl.Tx, tr.Tx, bl.Tx, tr.Tx]),
                           u.Quantity([bl.Ty, bl.Ty, tr.Ty, tr.Ty]),
                           frame=bl.frame.replicate_without_data())
        # On-disk coordinates are assumed to be on the solar surface, which is the
        # default sunpy assumption. Off-disk coordinates have no such natural
        # assumption, so they are placed on a spherical screen centered on the
        # original observer.
        with SphericalScreen(bl.frame.observer, only_off_disk=True):
            corners = corners.transform_to(earth_frame)

        if np.any(np.isnan(corners.Tx)) or np.any(np.isnan(corners.Ty)):
            raise ValueError("The requested cutout could not be transformed to a "
                             "Helioprojective frame with Earth as the observer. This "
                             "normally means part of the requested field of view is not "
                             "visible from Earth.")

        new_bl = SkyCoord(corners.Tx.min(), corners.Ty.min(), frame=earth_frame)
        new_tr = SkyCoord(corners.Tx.max(), corners.Ty.max(), frame=earth_frame)

        shift = np.sqrt((new_bl.Tx - bl.Tx)**2 + (new_bl.Ty - bl.Ty)**2).to('arcsec')
        warn_user(
            "The cutout coordinates have been transformed to a Helioprojective frame with "
            "Earth as the observer, which is used as an approximation of SDO's location. "
            f"The bottom-left corner moved by {shift:.3f}. On-disk coordinates were assumed "
            "to lie on the solar surface, and off-disk coordinates were assumed to lie on a "
            "spherical screen centered on the original observer. See the documentation for "
            "`sunpy.net.jsoc.Cutout` for details."
        )
        return new_bl, new_tr

    def collides(self, other):
        return isinstance(other, self.__class__)


walker = AttrWalker()


@walker.add_creator(AttrOr)
def _create1(wlk, query):
    qblocks = []
    for iattr in query.attrs:
        qblocks.extend(wlk.create(iattr))
    return qblocks


@walker.add_creator(AttrAnd, DataAttr)
def _create(wlk, query):
    map_ = {}
    wlk.apply(query, map_)
    return [map_]


@walker.add_applier(AttrAnd)
def _apply(wlk, query, imap):
    for iattr in query.attrs:
        wlk.apply(iattr, imap)


@walker.add_applier(SimpleAttr)
def _apply1(wlk, query, imap):
    imap[query.__class__.__name__.lower()] = query.value


@walker.add_applier(PrimeKey)
def _apply1(wlk, query, imap):
    key = 'primekey'
    if key in imap:
        imap[key][query.label] = query.value
    else:
        imap[key] = {query.label: query.value}


@walker.add_applier(Keyword)
def _apply1(wlk, query, imap):
    raise ValueError(f"Keyword '{query.value}' needs to have a comparison to a value.")


@walker.add_applier(KeywordComparison)
def _apply1(wlk, query, imap):
    key = 'keyword'
    if key in imap:
        imap[key][query.name] = {"operator": query.operator, "value": query.value}
    else:
        imap[key] = {f"{query.name}": {"operator": query.operator, "value": query.value}}


@walker.add_applier(Segment)
def _apply1(wlk, query, imap):
    key = 'segment'
    if key in imap:
        imap[key].append(query.value)
    else:
        imap[key] = [query.value]


@walker.add_applier(Cutout)
def _apply1(wlk, query, imap):
    imap[query.__class__.__name__.lower()] = query.value


@walker.add_applier(Time)
def _apply1(wlk, query, imap):
    imap['start_time'] = query.start
    imap['end_time'] = query.end


@walker.add_applier(Wavelength)
def _apply1(wlk, query, imap):
    if query.min != query.max:
        raise ValueError(
            "For JSOC queries Wavelength.min must equal Wavelength.max")

    imap[query.__class__.__name__.lower()] = query.min
