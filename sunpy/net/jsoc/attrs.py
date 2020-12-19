import astropy.units as u

from sunpy.coordinates.utils import get_rectangle_coordinates
from sunpy.net._attrs import Time, Wavelength
from sunpy.net.attr import AttrAnd, AttrOr, AttrWalker, DataAttr, SimpleAttr
from sunpy.util.decorators import deprecated

__all__ = ['Series', 'Protocol', 'Notify', 'Segment', 'Keys', 'PrimeKey',
           'Cutout']


# Define a custom __dir__ to restrict tab-completion to __all__
def __dir__():
    return __all__


class Series(SimpleAttr):
    """
    The JSOC Series to Download.

    This is the list of `Series <http://jsoc.stanford.edu/JsocSeries_DataProducts_map.html>`__.
    """


@deprecated(since="2.1", message="specify desired keywords as arguments to JSOCResponse.show()")
class Keys(SimpleAttr):
    """
    Keys choose which keywords to fetch while making a query request.
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
        Coordinate for the bottom left corner of the cutout.
    top_right : `~astropy.coordinates.SkyCoord`, optional
        Coordinate for the top right corner of the cutout. If this is
        not specified, both `width` and `height` must both be specified.
    width : `~astropy.units.Quantity`, optional
        Width of the cutout. If this parameter, along with `height`, is
        not specified, `top_right` must be specified.
    height : `~astropy.units.Quantity`, optional
        Height of the cutout. If this parameter, along with `width`, is
        not specified, `top_right` must be specified.
    tracking : `bool`, optional
        If True, the field of view follows the rotation of the Sun
    register : `bool`, optional
        If True, use sub-pixel registration when cropping to the target location.
    nan_off_limb : `bool`, optional
        If True, all off-limb pixels are set to NaN

    See Also
    --------
    sunpy.coordinate.utils.get_rectangle_coordinates
    """

    @u.quantity_input
    def __init__(self, bottom_left, top_right=None, width: u.arcsec = None,
                 height: u.arcsec = None, tracking=False, register=False,
                 nan_off_limb=False):
        super().__init__()
        bl, tr = get_rectangle_coordinates(bottom_left, top_right=top_right, width=width,
                                           height=height)
        self.value = {
            't_ref': bl.obstime.isot,
            # JSOC input is disable tracking so take the negative
            't': int(not tracking),
            'r': int(register),
            'c': int(nan_off_limb),
            'locunits': 'arcsec',
            'boxunits': 'arcsec',
            'x': ((bl.Tx + tr.Tx) / 2).to('arcsec').value,
            'y': ((bl.Ty + tr.Ty) / 2).to('arcsec').value,
            'width': (tr.Tx - bl.Tx).to('arcsec').value,
            'height': (tr.Ty - bl.Ty).to('arcsec').value,
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
