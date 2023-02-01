"""
===================================================================
Plotting of near real time (NRT) data by fetching it from json file
====================================================================

How to query and visualize nrt data(json format) using pandas and plot 
it using TimeSeries Library
"""



import pandas as pd 
from sunpy.time import parse_time
from sunpy import timeseries as ts 
from collections import OrderedDict
from astropy import units as u 

data = pd.read_json("https://services.swpc.noaa.gov/json/goes/primary/xrays-7-day.json")

############################################
#over here the json is transformed and stored as a pandas dataframe

data_short = data[data["energy"]=="0.05-0.4nm"]

##################################################################
#Filtering the data: A new dataframes are created by filtering the original data. 
# The data_short dataframe contains only the rows where the value of the "energy" column is "0.05-0.4nm"

data_long = data[data["energy"]=="0.1-0.8nm"]

##################################################################
# Again a new temporary dataframe is formed ,and the data_long contains only 
# the rows where the value of the "energy" column is "0.1-0.8nm".

time_array = parse_time(data_short["time_tag"])

##################################################
# The "time_tag" column from the data_short dataframe is passed 
# to the parse_time function to parse the strings into Time objects.

units = OrderedDict([("xrsa", u.W/u.m**2), ("xrsb", u.W/u.m**2)])

#####################################################
#The units variable is defined as an OrderedDict that maps the names of the two columns, 
# "xrsa" and "xrsb", to the corresponding physical units, u.W/u.m**2. 
# this plots data for energy range from data_short as xrsa and energy range from dta_long as xrsb

meta = OrderedDict({"instrument":"GOES X-ray sensor", "measurments": "primary", "type": "quicklook"}) 

#######################################################
# The meta variable is defined as an OrderedDict that contains information about the data, 
# such as the instrument used to obtain the data and the type of data.
#it can contain other types of metadata as well

goes_data = pd.DataFrame({"xrsa": data_short["flux"].values, "xrsb": data_long["flux"].values}, index=time_array.datetime)

###########################################################
# after filtering all the data the final data is merged into goes_data as final dataset
#Creating the time series data: A new Pandas dataframe, goes_data, is created 
# using the values of the "flux" columns from the data_short 
# and data_long dataframes, with the parsed time data as the index.

test_ts = ts.TimeSeries(goes_data, meta, units, source="xrs")

####################################################################
#Creating the time series object: The TimeSeries class from the sunpy.timeseries module is used to create a time series object, 
# test_ts, from the goes_data dataframe, the meta metadata, the units units, and a source string "xrs".

test_ts.peek()

############################################################
#Plotting the time series: The peek method of the test_ts object is called to plot the time series data.
