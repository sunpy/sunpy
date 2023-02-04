"""
==================================
How to create a `sunpy.timeseries.TimeSeries` from the GOES-XRS near real time (NRT). 
==================================

Currently, ``sunpy`` has no mechanism to directly download the GOES XRS Near Real Time (NRT) data.

We will showcase an example of how to download and load these files into a `sunpy.timeseries.TimeSeries`.
"""

import pandas as pd

from astropy import units as u

from sunpy import timeseries as ts
from sunpy.time import parse_time

###############################################################################
# We will start by getting reading the GOES JSON file using `pandas`.
# This allows us to download the file and load it straight into a `pandas.DataFrame`.
# This file updates every minute and contains only the last 7 days worth of data.

goes_json_data = pd.read_json("https://services.swpc.noaa.gov/json/goes/primary/xrays-7-day.json")

###############################################################################
# XRS collects data in two energy channels, "0.05-0.4nm" and "0.1-0.8nm".
# We separate these `short` and `long` wavelength readings into two arrays.

# This will get us the short wavelength data.
goes_short = data[data["energy"] == "0.05-0.4nm"]
# This will get us the long wavelength data.
goes_long = data[data["energy"] == "0.1-0.8nm"]

###############################################################################
# Similar to a `pandas.DataFrame`, `sunpy.timeseries.TimeSeries` requires
# a datetime index which we can get directly and transform into `astropy.time.Time`.

time_array = parse_time(data_short["time_tag"])

###############################################################################
# `sunpy.timeseries.TimeSeries` requires that there are units for data variables.
# To do this, we will create a dictionary that will map the names of the two columns,
# "xrsa" and "xrsb" (the channel names for GOES XRS), to their corresponding physical flux units, ``u.W/u.m**2``.

units = dict([("xrsa", u.W/u.m**2), ("xrsb", u.W/u.m**2)])

###############################################################################
# Typically `sunpy.timeseries.TimeSeries` will read metadata from the file,
# however, here we need to define our own metadata and we will keep it fairly simple.

meta = dict({"instrument": "GOES X-ray sensor", "measurements": "primary", "type": "quicklook"})

###############################################################################
# The meta variable is defined as an OrderedDict that contains information about the data,
# such as the instrument used to obtain the data and the type of data.
# it can contain other types of metadata as well

goes_data = pd.DataFrame({"xrsa": data_short["flux"].values, "xrsb": data_long["flux"].values}, index=time_array.datetime)

###############################################################################
# Now we will create a `sunpy.timeseries.TimeSeries` by passing in the data,
# the metadata and the units.

test_ts = ts.TimeSeries(goes_data, meta, units, source="xrs")

###############################################################################
# Finally, we can plot the timeseries.

test_ts.plot()

plt.show()
