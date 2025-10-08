import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.coordinates.frames import Helioprojective
from sunpy.coordinates.utils import get_rectangle_coordinates
from sunpy.map.maputils import coordinate_is_on_solar_disk
from sunpy.net._attrs import Time, Wavelength
from sunpy.net.attr import AttrAnd, AttrComparison, AttrOr, AttrWalker, DataAttr, SimpleAttr

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
    frame. The ``observer`` frame attribute of ``bottom_left`` is ignored, and
    instead the observer location is assumed to be SDO.

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

        center_x = (bl.Tx + tr.Tx) / 2
        center_y = (bl.Ty + tr.Ty) / 2
        center = SkyCoord(center_x, center_y, frame=bottom_left.frame)
        if tracking and not coordinate_is_on_solar_disk(center):
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
