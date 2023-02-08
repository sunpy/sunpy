"""
================================================
Getting the location of STEREO-A using solarmach
================================================

How to get the position of planetary bodies im the solar system using
`astropy's solar system ephemeris <http://docs.astropy.org/en/stable/coordinates/solarsystem.html#solar-system-ephemerides>`__ information and sunpy.
"""
import datetime

from solarmach import SolarMACH

##############################################################################
# These are necessary options

body_list = ['STEREO-A', 'Earth', 'Mars']
vsw_list = [400, 400, 400]
time = datetime.datetime.now()
date = str(time)

##############################################################################
# These are optional parameters

coord_sys = 'Stonyhurst'
reference_long = 273
reference_lat = 0
plot_spirals = True
plot_sun_body_line = True
long_offset = 270
reference_vsw = 400
return_plot_object = False
transparent = False
numbered_markers = True

##############################################################################
# initializing

sm = SolarMACH(date, body_list, vsw_list, reference_long, reference_lat, coord_sys)

##############################################################################
# Let's plot the results. Remember the Sun is at the center of this coordinate
# system.

sm.plot(
    plot_spirals=plot_spirals,
    plot_sun_body_line=plot_sun_body_line,
    reference_vsw=reference_vsw,
    transparent=transparent,
    numbered_markers=numbered_markers,
    long_offset=long_offset,
    return_plot_object=return_plot_object,

)
