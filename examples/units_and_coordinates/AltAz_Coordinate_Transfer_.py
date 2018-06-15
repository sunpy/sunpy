"""
============================================
Helioprojective/AltAz Coordinate Conversions
============================================

This is a very simple example to show how to convert
Helioprojective Coordiantes to Altitude/Azimuth Coordinates.
In this example, we use balloons observations to include
elevation as a factor.
"""
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
from astropy.time import Time
from sunpy.coordinates import frames, get_sunearth_distance
import astropy.units as u

######################################################################################
# We use `~astropy.coordinates.SkyCoord` to define the center of the Sun
obstime = "2013-09-21 16:00:00"
c = SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime=obstime, frame=frames.Helioprojective)

######################################################################################
# Now we establish our location on the Earth, in this case Fort Sumner, NM.
# We use the balloon's observational altitude as 'height'. Accuracy of 'height' is
# far less of a concern than Lon/Lat accuracy.
Fort_Sumner = EarthLocation(lat=34.4900*u.deg, lon=-104.221800*u.deg, height=40*u.km)

######################################################################################
# Now lets convert this to a local measurement of Altitude and Azimuth.

frame_altaz = AltAz(obstime=Time(obstime), location=Fort_Sumner)
sun_altaz = c.transform_to(frame_altaz)
print('Altitude is {0} and Azimuth is {1}'.format(sun_altaz.T.alt, sun_altaz.T.az))

######################################################################################
# Next let's check this calculation by converting it back to helioprojective.
# We should get our original input which was the center of the Sun.
# To go from Altitude/Azimuth to Helioprojective, you will need the distance to the Sun.
# solar distance. Define distance with SunPy's almanac.
distance = get_sunearth_distance(obstime)
b = SkyCoord(az=sun_altaz.T.az, alt=sun_altaz.T.alt, distance=distance, frame=frame_altaz)
sun_helio = b.transform_to(frames.Helioprojective)
print('The helioprojective point is {0}, {1}'.format(sun_helio.T.Tx, sun_helio.T.Ty))

######################################################################################
# The output is within a radius of 0.02 arcseccs.
