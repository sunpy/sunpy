"""
================================================
Getting the location of STEREO-A using solarmach
================================================

This example showcases how to get and plot the position of STEREO (Solar Terrestrial Relations Observatory) and other planetary bodies within the solar system using `solarmach <https://github.com/jgieseler/solarmach/>`__.
"""

import datetime

from solarmach import SolarMACH

##############################################################################
# These are the necessary arguments which are to be passed to ``solarMACH()``.
# You can make use of `~datetime.datetime.now` to pass current date and time.

body_list = ['STEREO-A', 'Earth', 'Mars']
vsw_list = [400, 400, 400]
time = datetime.datetime.now()
date = str(time)

##############################################################################
# If you want to make the plot look more clean and shart then you can pass
# these optional arguments as parameters to ``solarMACH()``.
# You can refer to `this <https://github.com/jgieseler/solarmach#readme>`__
# for docs.

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
# Now pass all the necessary and optional arguments as parameters to ``solarMACH()``.

sm = SolarMACH(date, body_list, vsw_list, reference_long, reference_lat, coord_sys)

##############################################################################
# Now plot the results. Remember the Sun is at the center of this coordinate
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
