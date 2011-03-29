# -*- coding: utf-8 -*-
#
# Author: Steven Christe <steven.d.christe@nasa.gov>
#
# <License info will go here...>
"""Sun-related parameters"""

import datetime
import math

def radius(t=None):
 """Returns the apparent radius of the Sun in arcsec.
 
    Based on IDL function get_sun.pro
    
    Notes:
    
    Keith (2011/03/29)
     If possible, it might be helpful to encode some of the magic
     numbers as constant, e.g. CONSTANT_NAME = 12345.6
    
     Also, it probably won't make too much of a difference, but
     you could either use numpy.pi, or defined pi once (with
     the same precision as numpy) at the top of the page.
     
     Feel free to delete these comments once you've read them :)

 """
 if t is None:
    t = datetime.datetime.now()
 
 #item passed in must be a datetime object. Need to check for this but how?
 #many of the following tasks should be broken out as other functions
 #t = datetime.datetime(2011,1,4,6,40,34)
 
 # Julian Centuries from 1900.0:
 jul = t - datetime.datetime(1900, 1, 1, 0, 0, 0)

 #jul.days needs to return fractional days but does not here
 #need to fix this
 t = jul.days / 36525.0

 #Julian date
 jd = jul.days + 2415020.5

 # Carrington Rotation Number:
 carr = (1. / 27.2753) * (jd - 2398167.0) + 1.0

 # Geometric Mean Longitude (deg):
 mnl = 279.696680 + 36000.76892 * t + 0.0003025 * t**2
 mnl = mnl % 360.0

 # Mean anomaly (deg):
 mna = 358.475830 + 35999.049750 * t - 0.0001500 * t**2 - 0.00000330 * t**3
 mna = mna % 360.0

 # Eccentricity of orbit:
 e = 0.016751040 - 0.00004180 * t - 0.0000001260 * t**2

 # Sun's equation of center (deg):
 c = ((1.9194600 - 0.0047890 * t - 0.0000140 * t**2) 
   * math.sin(mna * 3.1415 / 180.0) + (0.0200940 - 0.0001000 * t) 
   * math.sin(2 * mna * 3.1415 / 180.0) + 0.0002930 
   * math.sin(3 * mna * 3.1415 / 180.0))

 # Sun's true geometric longitude (deg)
 #   (Refered to the mean equinox of date.  Question: Should the higher
 #    accuracy terms from which app_long is derived be added to true_long?)
 true_long = (mnl + c) % 360.0

 # Sun's true anomaly (deg):
 ta = (mna + c) % 360.0

 # Sun's radius vector (AU).  There are a set of higher accuracy
 #   terms not included here.  The values calculated here agree with
 #   the example in the book:
 dist = 1.00000020 * (1.0 - e**2) / (1.0 + e * math.cos(ta * 3.1415 / 180.0))

 # Semidiameter (arc sec):
 semidiameter_arsec = 959.63 / dist

 return semidiameter_arsec

def mass(units = 'g'):
 """Returns the mass of the Sun."""
 #would be nice to have some sort of units module to handle conversions
 m = 1.9891e33
 
 if units is 'g':
    return m
 elif units is 'kg':
    return m / 1000.0
 if units is 'earth':
    return m / 5.9742e24

