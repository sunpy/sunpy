"""
===================================================================
Plotting of near real time (NRT) data by fetching it from json file
====================================================================

How to query and visualize nrt data(json format) using pandas and plot
it using TimeSeries Library
"""


from collections import OrderedDict

import pandas as pd

from astropy import units as u

from sunpy import timeseries as ts
from sunpy.time import parse_time

data = pd.read_json("https://services.swpc.noaa.gov/json/goes/primary/xrays-7-day.json")

############################################
# over here the json is transformed and stored as a pandas dataframe

data_short = data[data["energy"] == "0.05-0.4nm"]

##################################################################
# Filtering the data: A new dataframes are created by filtering the original data.
# The data_short dataframe contains only the rows where the value of the "energy" column is "0.05-0.4nm"

data_long = data[data["energy"] == "0.1-0.8nm"]

##################################################################
# Again a new temporary dataframe is formed ,and the data_long contains only
# the rows where the value of the "energy" column is "0.1-0.8nm".

time_array = parse_time(data_short["time_tag"])

##################################################
# The "time_tag" column from the data_short dataframe is passed
# to the parse_time function to parse the strings into Time objects.

units = OrderedDict([("xrsa", u.W/u.m**2), ("xrsb", u.W/u.m**2)])

#####################################################
# The units variable is defined as an OrderedDict that maps the names of the two columns,
# "xrsa" and "xrsb", to the corresponding physical units, u.W/u.m**2.
# this plots data for energy range from data_short as xrsa and energy range from dta_long as xrsb

meta = OrderedDict({"instrument": "GOES X-ray sensor", "measurments": "primary", "type": "quicklook"})

#######################################################
# The meta variable is defined as an OrderedDict that contains information about the data,
# such as the instrument used to obtain the data and the type of data.
# it can contain other types of metadata as well

goes_data = pd.DataFrame({"xrsa": data_short["flux"].values, "xrsb": data_long["flux"].values}, index=time_array.datetime)

###########################################################
# Now we will create a `sunpy.timeseries.TimeSeries` by passing in the data, the metadata and the units.

test_ts = ts.TimeSeries(goes_data, meta, units, source="xrs")

####################################################################
# Now we have created the timeseries, we will finally plot it.

test_ts.plot()

plt.show()
