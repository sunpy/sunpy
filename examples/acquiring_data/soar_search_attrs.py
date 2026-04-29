"""
===========================
Available search attributes
===========================

This example demonstrates the available search attributes
for SOAR currently supported by ``sunpy-soar``.
"""

import sunpy.net.attrs as a

#####################################################
# The easiest way to access search attributes is using
# the attribute registry provided by `sunpy.net.attrs`.
#
# When constructing a search for SOAR ``a.Time`` must be provided.
# Other search attributes can be used too - ``sunpy-soar`` recognizes the following:
#
# The third ``near`` argument to ``a.Time`` is not currently supported.
# You will have to manually filter the results if you want to find the one closest to a given time.
#
# For instrument the following are supported:
#
# - "EPD": "Energetic Particle Detector"
# - "EUI": "Extreme UV Imager"
# - "MAG": "Magnetometer"
# - "METIS": "Metis: Visible light and ultraviolet coronagraph"
# - "PHI": "Polarimetric and Helioseismic Imager"
# - "RPW": "Radio and Plasma Wave instrument"
# - "SOLOHI": "Solar Orbiter Heliospheric Imager"
# - "SPICE": "SPectral Investigation of the Coronal Environment"
# - "STIX": "Spectrometer Telescope for Imaging X-rays"
# - "SWA": "Solar Wind Analyser"
#
# For level the following are supported:
# L0, L1, L2, L3, LL01, LL02, LL03
#
# For product:
a.soar.Product
#####################################################
# For specific instrument detectors or sensors, see ``a.Detector``.
# However, some SOAR products require the use of ``a.soar.Sensor`` attribute instead:
a.soar.Sensor
