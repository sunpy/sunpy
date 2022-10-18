************************
Querying Helioviewer.org
************************

The Helioviewer Project now maintains a Python wrapper around their API.
It is called `hvpy <https://hvpy.readthedocs.io/en/latest/>`__ and their documentation has several examples of how to use it.

Migration guide
===============

If you are using ``HelioviewerClient`` this guide shows to switch to using ``hvpy``.

``HelioviewerClient().data_sources`` is replaced by `hvpy.DataSource`, an `~enum.Enum` that lists every known data source for Helioviewer.org.

.. code-block:: python

   >>> import hvpy
   >>> print(hvpy.DataSource.AIA_131)
   DataSource.AIA_131

It also supports tab-complete, to find the data source you want.

``HelioviewerClient().get_data_sources`` is replaced by :func:`hvpy.getDataSources`.

.. code-block:: python

   >>> import hvpy
   >>> hvpy.getDataSources()  # doctest: +REMOTE_DATA
   {'SDO': {'HMI': {'continuum': {'sourceId': 18,
   ...

``HelioviewerClient().get_closest_image`` is replaced by :func:`hvpy.getClosestImage`.

.. code-block:: python

   >>> from datetime import datetime
   >>> import hvpy
   >>> hvpy.getClosestImage(date=datetime.now(), sourceId=hvpy.DataSource.AIA_171)  # doctest: +REMOTE_DATA
   {'id': '13326...',...

``HelioviewerClient().download_jp2`` is replaced by :func:`hvpy.getJP2Image`, but you will have to wrap the call using :func:`hvpy.save_file` to save the data to disk.

.. code-block:: python

   >>> from datetime import datetime
   >>> import hvpy
   >>> filepath = hvpy.save_file(hvpy.getJP2Image(date=datetime.today(), sourceId=hvpy.DataSource.AIA_171), filename="~/example.jpeg")  # doctest: +REMOTE_DATA
   >>> filepath.unlink()

``HelioviewerClient().get_jp2_header`` is replaced by :func:`hvpy.getJP2Header`.
However you will need to make a call to :func:`hvpy.getClosestImage` to get the ID required.
Furthermore, the header is returned as a XML string, which you will need to parse.

.. code-block:: python

   >>> from datetime import datetime
   >>> import hvpy
   >>> metadata = hvpy.getClosestImage(date=datetime.today(), sourceId=hvpy.DataSource.AIA_171)  # doctest: +REMOTE_DATA
   >>> hvpy.getJP2Header(metadata['id'])  # doctest: +REMOTE_DATA
   '<?xml version="1.0" encoding="utf-8"?>...

``HelioviewerClient().download_png`` is replaced by `hvpy.createScreenshot`, it takes the same arguments as the old method expect for ``progress`` and ``directory`` which do not exist and adds ``filename`` so one is able to save the file, otherwise it will save it in the current working directory.

.. code-block:: python

   >>> from datetime import datetime
   >>> import hvpy
   >>> screenshot_location = hvpy.createScreenshot(
   ...     date=datetime.today(),
   ...     layers=hvpy.create_layers([(hvpy.DataSource.AIA_171, 100)]),
   ...     events=hvpy.create_events([hvpy.EventType.ACTIVE_REGION]),
   ...     eventLabels=True,
   ...     imageScale=1,
   ...     x0=0,
   ...     y0=0,
   ...     width=100,
   ...     height=100,
   ...     filename="my_screenshot",
   ... )  # doctest: +REMOTE_DATA
   >>> screenshot_location.unlink()  # doctest: +REMOTE_DATA

`The documentation for hvpy has more examples of how to use it and examples for each function <https://hvpy.readthedocs.io/en/latest/index.html>`__.

If you encounter a problem with the new API, please open an issue on `the hvpy issue tracker <https://github.com/Helioviewer-Project/python-api/issues>`__.
