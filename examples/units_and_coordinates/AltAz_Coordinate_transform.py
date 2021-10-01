"""
=======================================================
Converting between Helioprojective and AltAz Coordinate
=======================================================

How to find the Sun in the sky as viewed from a particular location.
"""
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time

from sunpy.coordinates import frames, sun

######################################################################################
# We use `~astropy.coordinates.SkyCoord` to define the center of the Sun.

obstime = "2013-09-21 16:00:00"
c = SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime=obstime,
             observer="earth", frame=frames.Helioprojective)

######################################################################################
# Now we establish our location on the Earth, in this case let's consider a high altitude
# balloon launched from Fort Sumner, NM.

Fort_Sumner = EarthLocation(lat=34.4900*u.deg, lon=-104.221800*u.deg, height=40*u.km)

######################################################################################
# Now lets convert this to a local measurement of Altitude and Azimuth.

frame_altaz = AltAz(obstime=Time(obstime), location=Fort_Sumner)
sun_altaz = c.transform_to(frame_altaz)
print(f'Altitude is {sun_altaz.T.alt} and Azimuth is {sun_altaz.T.az}')

######################################################################################
# Next let's check this calculation by converting it back to helioprojective.
# We should get our original input which was the center of the Sun.
# To go from Altitude/Azimuth to Helioprojective, you will need the distance to the Sun.
# solar distance. Define distance with sunpy's almanac.

distance = sun.earth_distance(obstime)
b = SkyCoord(az=sun_altaz.T.az, alt=sun_altaz.T.alt, distance=distance, frame=frame_altaz)
sun_helio = b.transform_to(frames.Helioprojective(observer="earth"))
print(f'The helioprojective point is {sun_helio.T.Tx}, {sun_helio.T.Ty}')
