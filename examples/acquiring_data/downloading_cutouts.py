from aiapy.util import sdo_location

import astropy.time
import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.coordinates import Helioprojective
from sunpy.net import attrs, jsoc

tstart = astropy.time.Time('2015-10-17T04:33:30.000')
duration = 12 * u.h

obs = sdo_location(tstart)
frame = Helioprojective(observer=obs, obstime=obs.obstime)

width = 345.6 * u.arcsec
height = 345.6 * u.arcsec
Txc = -517.2 * u.arcsec
Tyc = -246 * u.arcsec

bottom_left = SkyCoord(Tx=Txc-width/2, Ty=Tyc-height/2, frame=frame)

c = jsoc.JSOCClient()

q = c.search(
    attrs.Time(tstart, tstart + duration),
    attrs.Wavelength(171*u.angstrom),
    attrs.Sample(12*u.min),
    jsoc.attrs.Series('aia.lev1_euv_12s'),
    jsoc.attrs.Notify('will.t.barnes@gmail.com'),
    jsoc.attrs.Segment('image'),
    jsoc.attrs.Cutout(bottom_left, width=width, height=height, tracking=True),
)
resp = c.request_data(q, method='url', protocol='fits')
while not resp.has_succeeded():
    continue
req = c.get_request(
    resp,
    path='/Users/willbarnes/Desktop/cutout-request-tests', max_conn=2,
)
