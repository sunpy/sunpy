.. _sunpy-soar-dev-guide-working:


How data is retrieved from the SOAR
===================================

To retrieve data from SOAR, you first use the `sunpy.net.Fido` object from `sunpy` and specify the desired attributes using `sunpy.net.attrs`.
These attributes define the criteria for the data you want to retrieve, such as the time range, instrument, or wavelength.

Here is an example of how to specify the time range for the data you want to retrieve:

.. code-block:: python

    import sunpy.net.attrs as a
    from sunpy.net import Fido

    instrument = a.Instrument("EUI")
    time = a.Time("2021-02-01", "2021-02-02")
    level = a.Level(1)
    product = a.soar.Product("EUI-FSI174-IMAGE")

    result = Fido.search(instrument & time & level & product)


``sunpy-soar`` constructs the query based on the specified criteria, then generates a URL to interact with the SOAR API.
This is done by using the Table Access Protocol (TAP), a widely adopted standard in the astronomical community for accessing large datasets.
The queries are formulated in Astronomical Data Query Language (`ADQL <https://www.ivoa.net/documents/ADQL/>`__), which allows for flexible querying of astronomical data.
The results are returned in the form of an Astropy table, providing a structured and efficient format for further analysis and visualization within the Python environment.

A generated query looks like:

.. code-block:: SQL

    SELECT * FROM v_sc_data_item WHERE instrument='EPD' AND begin_time>='2021-02-01 00:00:00' AND begin_time<='2021-02-02 00:00:00' AND level='L1' AND descriptor='epd-epthet2-nom-close'

Or with a JOIN

.. code-block:: SQL

    SELECT h1.instrument, h1.descriptor, h1.level, h1.begin_time, h1.end_time, h1.data_item_id, h1.filesize, h1.filename, h1.soop_name, h2.detector, h2.wavelength, h2.dimension_index
    FROM v_sc_data_item AS h1 JOIN v_eui_sc_fits AS h2 USING (data_item_oid) WHERE h1.instrument='EUI' AND h1.begin_time>='2021-02-01 00:00:00' AND h1.begin_time<='2021-02-02 00:00:00' AND
    h2.dimension_index='1' AND h1.level='L1' AND h1.descriptor='eui-fsi174-image'

The URL is generated with the query formed based on the parameters, then Fido is used to search and download the data.
