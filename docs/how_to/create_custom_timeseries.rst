.. _sunpy-how-to-custom-timeseries:

************************
Create Custom TimeSeries
************************

Sometimes you will have data that you want to transform into a TimeSeries.
You can use the factory to create a `~sunpy.timeseries.GenericTimeSeries` from a variety of data sources currently including `pandas.DataFrame` and `astropy.table.Table`.

Creating a TimeSeries from a Pandas DataFrame
=============================================

A TimeSeries object must be supplied with some data when it is created.
The data can either be in your current Python session, in a local file, or in a remote file.
Let's create some data and pass it into a TimeSeries object:

.. code-block:: python

    >>> import numpy as np

    >>> intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

This creates a basic numpy array of values representing a sine wave.
We can use this array along with a suitable time storing object (such as Astropy `~astropy.time` or a list of `datetime` objects) to make a Pandas `~pandas.DataFrame`.
A suitable list of times must contain the same number of values as the data, this can be created using:

.. code-block:: python

    >>> import datetime

    >>> base = datetime.datetime.today()
    >>> times = [base - datetime.timedelta(minutes=x) for x in range(24*60, 0, -1)]

The Pandas `~pandas.DataFrame` will use the dates list as the index:

.. code-block:: python

    >>> from pandas import DataFrame

    >>> data = DataFrame(intensity, index=times, columns=['intensity'])

This `~pandas.DataFrame` can then be used to construct a TimeSeries:

.. code-block:: python

    >>> import astropy.units as u

    >>> import sunpy.timeseries as ts

    >>> header = {'key': 'value'}
    >>> units = {'intensity': u.W/u.m**2}
    >>> ts_custom = ts.TimeSeries(data, header, units)

Creating Custom TimeSeries from an Astropy Table
================================================

A Pandas `~pandas.DataFrame` is the underlying object used to store the data within a TimeSeries, so the above example is the most lightweight to create a custom TimeSeries, but being scientific data it will often be more convenient to use an Astropy `~astropy.table.Table` to create a TimeSeries.
An advantage of this method is it allows you to include metadata and Astropy `~astropy.units.quantity.Quantity` values, which are both supported in tables, without additional arguments.
For example:

.. code-block:: python

    >>> from astropy.time import Time
    >>> from astropy.table import Table

    >>> base = datetime.datetime.today()
    >>> times = [base - datetime.timedelta(minutes=x) for x in range(24*60, 0, -1)]
    >>> intensity = u.Quantity(np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60)))), u.W/u.m**2)
    >>> tbl_meta = {'t_key':'t_value'}
    >>> table = Table([times, intensity], names=['time', 'intensity'], meta=tbl_meta)
    >>> table.add_index('time')
    >>> ts_table = ts.TimeSeries(table)

Note that due to the properties of the `~astropy.time.Time` object, this will be a mixin column which since it is a single object, limits the versatility of the `~astropy.table.Table` a little.
For more on mixin columns see the `Astropy docs <https://docs.astropy.org/en/stable/table/mixin_columns.html>`__.
The units will be taken from the table quantities for each column, the metadata will simply be the ``table.meta`` dictionary.
You can also explicitly add metadata and units, these will be added to the relevant dictionaries using the dictionary update method, with the explicit user-given values taking precedence:

.. code-block:: python

    >>> from collections import OrderedDict

    >>> from sunpy.util.metadata import MetaDict

    >>> meta = MetaDict({'key':'value'})
    >>> units = OrderedDict([('intensity', u.W/u.m**2)])
    >>> ts_table = ts.TimeSeries(table, meta, units)
