"""
================================================
Getting the location of STEREO-A using solarmach
================================================

This example showcases how to get and plot the position of STEREO (Solar Terrestrial Relations Observatory) and other planetary bodies within the solar system using `solarmach <https://github.com/jgieseler/solarmach/>`__.
"""

import datetime

from solarmach import SolarMACH

##############################################################################
# Since we will use ``solarmach``, we will need to define the required paramters:
# You can find an example of more complex use on the `GitHub readme <https://github.com/jgieseler/solarmach#usage>`__.
# For now we will define the bodies we want, the position-sensitive solar wind speed and the time for the image.

bodies = ['STEREO-A', 'Earth', 'Mars']
wind_speeds = [400, 400, 400]
today = str(datetime.datetime.now())

##############################################################################
##############################################################################
# ``solarmach.SolarMACH`` has many optional parameters and they are detailed in `GitHub readme <https://github.com/jgieseler/solarmach#usage>`__.
# For this example, the defaults fit exactly what we are after:

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
