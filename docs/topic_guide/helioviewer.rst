.. _sunpy-topic-guide-using-helioviewer:

****************************************
Switching from HelioviewerClient to hvpy
****************************************

In older versions of sunpy it was possible to query Helioviewer using the ``HelioviewerClient``.
The Helioviewer Project now maintains a Python wrapper around their API that supersedes this.
It is called `hvpy <https://hvpy.readthedocs.io/en/latest/>`__ and their documentation has several examples of how to use it.
This page explains how to migrate from ``HelioviewerClient`` to ``hvpy``.

``HelioviewerClient().data_sources`` is replaced by `hvpy.DataSource`, an `~enum.Enum` that lists every known data source for Helioviewer.org.

.. code-block:: python

   >>> import hvpy

   >>> print(hvpy.DataSource.AIA_131)
   DataSource.AIA_131

It also supports tab-complete, to find the data source you want.

``HelioviewerClient().get_data_sources`` is replaced by :func:`hvpy.getDataSources`.

.. code-block:: python

   >>> hvpy.getDataSources()  # doctest: +REMOTE_DATA
   {...'SDO': {...'HMI': {'continuum': {'sourceId': 18,
   ...

``HelioviewerClient().get_closest_image`` is replaced by :func:`hvpy.getClosestImage`.

.. code-block:: python

   >>> from datetime import datetime

   >>> hvpy.getClosestImage(date=datetime(2022, 1, 1), sourceId=hvpy.DataSource.AIA_171)  # doctest: +REMOTE_DATA
    {'id': '79024526', 'date': '2022-01-01 00:04:57', 'name': 'AIA 171', 'scale': 0.5899466652089547, 'scaleCorrection': 1.0170411248743723, 'width': 4096, 'height': 4096, 'refPixelX': 2048.5, 'refPixelY': 2048.5, 'rotation': 0, 'rsun': 1626.6638, 'dsun': 147091270000, 'sunCenterOffsetParams': [], 'layeringOrder': 1}

``HelioviewerClient().download_jp2`` is replaced by :func:`hvpy.getJP2Image`, but you will have to wrap the call using :func:`hvpy.save_file` to save the data to disk.

.. code-block:: python

   >>> filepath = hvpy.save_file(hvpy.getJP2Image(date=datetime.today(), sourceId=hvpy.DataSource.AIA_171), filename="~/example.jpeg")  # doctest: +REMOTE_DATA
   >>> filepath.unlink()  # doctest: +REMOTE_DATA

``HelioviewerClient().get_jp2_header`` is replaced by :func:`hvpy.getJP2Header`.
However you will need to make a call to :func:`hvpy.getClosestImage` to get the ID required.
Furthermore, the header is returned as a XML string, which you will need to parse.

.. code-block:: python

   >>> metadata = hvpy.getClosestImage(date=datetime.today(), sourceId=hvpy.DataSource.AIA_171)  # doctest: +REMOTE_DATA
   >>> hvpy.getJP2Header(metadata['id'])  # doctest: +REMOTE_DATA
   '<?xml version="1.0" encoding="utf-8"?>...

``HelioviewerClient().download_png`` is replaced by `hvpy.createScreenshot`, it takes the same arguments as the old method expect for ``progress`` and ``directory`` which do not exist and adds ``filename`` so one is able to save the file, otherwise it will save it in the current working directory.

.. code-block:: python

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
