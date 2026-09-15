"""
========================
Getting data from CDAWeb
========================

How to download data from the Coordinated Data Analysis Web (CDAWeb).

CDAWeb stores data from from current and past space physics missions, and is
full of heliospheric insitu datasets.
"""
# sphinx_gallery_tags = ["Acquiring Data", "CDAWeb", "Solar Orbiter"]

import itables
from IPython.display import HTML

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net.attr import _create_table
from sunpy.timeseries import TimeSeries

###############################################################################
# `sunpy.net.Fido` is the primary interface to search for and download data and
# will automatically search CDAWeb when the ``cdaweb.Dataset`` attribute is provided to
# the search. To lookup the different dataset IDs available, you can use the
# form at https://cdaweb.gsfc.nasa.gov/index.html/
#
# There are thousands of dataset names, so it's usually much quicker to search
# through them from Python rather than the form above. If you're working in a
# Jupyter notebook and have the optional ``itables`` package installed, you can
# pull up every dataset name in an interactive, filterable table with
# ``a.cdaweb.Dataset.show_in_notebook()``.
#
# The same table is shown below, and you can type into the search box to
# filter it down to the dataset you're after. ``show_in_notebook`` only
# displays a table when it's called inside a live notebook, so here we build
# the same table by hand from the underlying attr registry, and bump
# ``maxBytes`` so none of the ~3000 rows get dropped.
itables.options.maxBytes = "1MB"
HTML(itables.to_html_datatable(_create_table(a.cdaweb.Dataset).to_pandas()))

###############################################################################
# Once you've found the dataset you want, you can pass its name straight to
# `~sunpy.net.attrs.cdaweb.Dataset`.
trange = a.Time('2021/07/01', '2021/07/08')
dataset = a.cdaweb.Dataset('SOLO_L2_MAG-RTN-NORMAL-1-MINUTE')
result = Fido.search(trange, dataset)

###############################################################################
# Let's inspect the results. We can see that there's seven files, one for each
# day within the query.
print(result)

###############################################################################
# Let's download the first two files
downloaded_files = Fido.fetch(result[0, 0:2])
print(downloaded_files)

###############################################################################
# Finally we can load and take a look at the data using
# `~sunpy.timeseries.TimeSeries` This requires an installation of the cdflib
# Python library to read the CDF file.
solo_mag = TimeSeries(downloaded_files, concatenate=True)
print(solo_mag.columns)
solo_mag.peek(columns=['B_RTN_0', 'B_RTN_1', 'B_RTN_2'])
