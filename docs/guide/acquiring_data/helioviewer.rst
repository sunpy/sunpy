-----------------------------------
Querying Helioviewer.org with SunPy
-----------------------------------
SunPy can be used to make several basic requests using the The `Helioviewer.org API <https://api.helioviewer.org/docs/v2/>`_
including generating a PNG and downloading a `JPEG 2000 <http://wiki.helioviewer.org/wiki/JPEG_2000>`_
image and loading it into a SunPy Map.

The SunPy Helioviewer client requires installing two other pieces of software.
The first OpenJPEG is an open source library for reading and writing JPEG2000
files.  To install OpenJPEG, please follow the instructions at `the OpenJPEG
homepage <http://www.openjpeg.org>`_.

The other package you will need is `Glymur
<https://pypi.python.org/pypi/Glymur/>`_.  Glymur is an interface
between Python and the OpenJPEG libraries.  Please follow the
instructions `here <https://glymur.readthedocs.io/en/latest/>`_ to
install Glymur on your system.

To interact with the Helioviewer API, users first create a "HelioviewerClient"
instance. The client instance can then be used to make various queries against
the API using the same parameters one would use when making a web request. Note that
the HelioviewerClient does not currently offer full access to the HelioViewer API. We provide functions:

1. That download a JP2 image for the specified datasource that is the closest match in time to the `date` requested.
2. That returns a hierarchial list of the available datasources.
3. To find the image data that is closest to the requested date/time. This function returns the associated metadata from the helioviewer database and the XML header of the JPEG2000 image file and,
4. To generate custom screenshots.

Nearly all requests require the user to specify the data they are interested in
and this can be done using the following methods:

1. Use the datasource_info dictionary to get the list of data source and their corresponding source id
numbers which can be used to refer to a particular dataset or,
2. Specify the four components of a Helioviewer.org data source or layer: *observatory*, *instrument*, 
*detector* and *measurement*. The default value is *None* for all the components and should be changed
to the correct value for the source you are interested in.

Let us begin by retrieving the available list of sources that Helioviewer supports by using `~sunpy.net.helioviewer.HelioviewerClient.datasource_info`:: 

    >>> from sunpy.net import helioviewer 
    >>> hv = helioviewer.HelioviewerClient()  # doctest: +REMOTE_DATA
    >>> datasource = hv.datasource_info  # doctest: +REMOTE_DATA
    >>> for sorted_ids in sorted(datasource.items(), key=lambda x: x[1]):  # doctest: +REMOTE_DATA
    ...     print(sorted_ids)  # doctest: +REMOTE_DATA
    (('SOHO', 'EIT', None, '171'), 0)
    (('SOHO', 'EIT', None, '195'), 1)
    (('SOHO', 'EIT', None, '284'), 2)
    (('SOHO', 'EIT', None, '304'), 3)
    (('SOHO', 'LASCO', 'C2', 'white-light'), 4)
    (('SOHO', 'LASCO', 'C3', 'white-light'), 5)
    (('SOHO', 'MDI', None, 'magnetogram'), 6)
    (('SOHO', 'MDI', None, 'continuum'), 7)
    (('SDO', 'AIA', None, '94'), 8)
    (('SDO', 'AIA', None, '131'), 9)
    (('SDO', 'AIA', None, '171'), 10)
    (('SDO', 'AIA', None, '193'), 11)
    (('SDO', 'AIA', None, '211'), 12)
    (('SDO', 'AIA', None, '304'), 13)
    (('SDO', 'AIA', None, '335'), 14)
    (('SDO', 'AIA', None, '1600'), 15)
    (('SDO', 'AIA', None, '1700'), 16)
    (('SDO', 'AIA', None, '4500'), 17)
    (('SDO', 'HMI', None, 'continuum'), 18)
    (('SDO', 'HMI', None, 'magnetogram'), 19)
    (('STEREO_A', 'SECCHI', 'EUVI', '171'), 20)
    (('STEREO_A', 'SECCHI', 'EUVI', '195'), 21)
    (('STEREO_A', 'SECCHI', 'EUVI', '284'), 22)
    (('STEREO_A', 'SECCHI', 'EUVI', '304'), 23)
    (('STEREO_B', 'SECCHI', 'EUVI', '171'), 24)
    (('STEREO_B', 'SECCHI', 'EUVI', '195'), 25)
    (('STEREO_B', 'SECCHI', 'EUVI', '284'), 26)
    (('STEREO_B', 'SECCHI', 'EUVI', '304'), 27)
    (('STEREO_A', 'SECCHI', 'COR1', 'white-light'), 28)
    (('STEREO_A', 'SECCHI', 'COR2', 'white-light'), 29)
    (('STEREO_B', 'SECCHI', 'COR1', 'white-light'), 30)
    (('STEREO_B', 'SECCHI', 'COR2', 'white-light'), 31)
    (('PROBA2', 'SWAP', None, '174'), 32)
    (('Yohkoh', 'SXT', None, 'AlMgMn'), 33)
    (('Yohkoh', 'SXT', None, 'thin-Al'), 34)
    (('Yohkoh', 'SXT', None, 'white-light'), 35)
    (('Hinode', 'XRT', 'Al_med', 'Al_thick'), 39)
    (('Hinode', 'XRT', 'Al_med', 'Be_thick'), 40)
    (('Hinode', 'XRT', 'Al_med', 'Open'), 42)
    (('Hinode', 'XRT', 'Al_med', 'Ti_poly'), 43)
    (('Hinode', 'XRT', 'Al_poly', 'Al_mesh'), 44)
    (('Hinode', 'XRT', 'Al_poly', 'Al_thick'), 45)
    (('Hinode', 'XRT', 'Al_poly', 'Be_thick'), 46)
    (('Hinode', 'XRT', 'Al_poly', 'Open'), 48)
    (('Hinode', 'XRT', 'Al_poly', 'Ti_poly'), 49)
    (('Hinode', 'XRT', 'Be_med', 'Open'), 54)
    (('Hinode', 'XRT', 'Be_thin', 'Open'), 60)
    (('Hinode', 'XRT', 'C_poly', 'Al_mesh'), 62)
    (('Hinode', 'XRT', 'C_poly', 'Al_thick'), 63)
    (('Hinode', 'XRT', 'C_poly', 'Open'), 66)
    (('Hinode', 'XRT', 'C_poly', 'Ti_poly'), 67)
    (('Hinode', 'XRT', 'Open', 'Al_mesh'), 69)
    (('Hinode', 'XRT', 'Open', 'Al_thick'), 70)
    (('Hinode', 'XRT', 'Open', 'Be_thick'), 71)
    (('Hinode', 'XRT', 'Open', 'Ti_poly'), 74)
    (('TRACE', None, None, '171'), 75)
    (('TRACE', None, None, '195'), 76)
    (('TRACE', None, None, '284'), 77)
    (('TRACE', None, None, '1216'), 78)
    (('TRACE', None, None, '1550'), 79)
    (('TRACE', None, None, '1600'), 80)
    (('TRACE', None, None, '1700'), 81)
    (('TRACE', None, None, 'white-light'), 82)
    (('Hinode', 'XRT', 'Any', 'Any'), 10001)
    (('Hinode', 'XRT', 'Any', 'Al_mesh'), 10002)
    (('Hinode', 'XRT', 'Any', 'Al_thick'), 10003)
    (('Hinode', 'XRT', 'Any', 'Be_thick'), 10004)
    (('Hinode', 'XRT', 'Any', 'Gband'), 10005)
    (('Hinode', 'XRT', 'Any', 'Open'), 10006)
    (('Hinode', 'XRT', 'Any', 'Ti_poly'), 10007)
    (('Hinode', 'XRT', 'Al_med', 'Any'), 10008)
    (('Hinode', 'XRT', 'Al_poly', 'Any'), 10009)
    (('Hinode', 'XRT', 'Be_med', 'Any'), 10010)
    (('Hinode', 'XRT', 'Be_thin', 'Any'), 10011)
    (('Hinode', 'XRT', 'C_poly', 'Any'), 10012)
    (('Hinode', 'XRT', 'Open', 'Any'), 10013)


Helioviewer provides JP2 images from a range of sources. New sources of JP2 images are being added every few months.

Suppose we next want to download a PNG image of the latest
AIA 304 image available on Helioviewer.org. We could use the explicit
approach as shown in the following example.::

   >>> from sunpy.net.helioviewer import HelioviewerClient
   >>> import matplotlib.pyplot as plt
   >>> from matplotlib.image import imread
   >>> hv = HelioviewerClient()  # doctest: +REMOTE_DATA
   >>> file = hv.download_png('2099/01/01', 4.8, "[SDO,AIA,AIA,304,1,100]", x0=0, y0=0, width=512, height=512)  # doctest: +REMOTE_DATA
   >>> im = imread(file)  # doctest: +REMOTE_DATA
   >>> plt.imshow(im)  # doctest: +SKIP
   >>> plt.axis('off')  # doctest: +SKIP
   >>> plt.show()  # doctest: +SKIP


.. image:: helioviewer-1.png


Where 4.8 refers to the image resolution in arcseconds per pixel (larger values
mean lower resolution), the "1" and "100" in the layer string refer to the
visibility (visible/hidden) and opacity, x0 and y0 are the center points about
which to focus and the width and height are the pixel values for the image
dimensions.

Note that the filename of the returned file has the date and time of
the request, not of any of the times shown in the image itself.  This
is not a bug.  Helioviewer.org finds images *closest to the requested
time*.  Since the user may ask for images from multiple sources, and
each of them may have a different observation time, the problem
becomes which time is the most appropriate to associate with the
resultant image.  Helioviewer.org doesn't choose between the images
times, but instead uses the request time to construct the image
filename.  This means that the image file names for request times in
the future (like in this example) can look a little unusual compared to
the times in the image.

If we find that the source id for AIA 304 is is 13, we could make the same
request using: ::

    hv.download_png('2099/01/01', 4.8, "[13,1,100]", x0=0, y0=0, width=512, height=512)

Now suppose we wanted to create a composite PNG image using data from two
different AIA wavelengths and LASCO C2 coronagraph data. The layer string is
extended to include the additional data sources, and opacity is throttled
down for the second AIA layer so that it does not completely block out the
lower layer.::

   >>> from sunpy.net.helioviewer import HelioviewerClient
   >>> import matplotlib.pyplot as plt
   >>> from matplotlib.image import imread
   >>> hv = HelioviewerClient()  # doctest: +REMOTE_DATA
   >>> file = hv.download_png('2099/01/01', 6, "[SDO,AIA,AIA,304,1,100],[SDO,AIA,AIA,193,1,50],[SOHO,LASCO,C2,white-light,1,100]", x0=0, y0=0, width=768, height=768)  # doctest: +REMOTE_DATA
   >>> im = imread(file)  # doctest: +REMOTE_DATA
   >>> plt.imshow(im)  # doctest: +SKIP
   >>> plt.axis('off')  # doctest: +SKIP
   >>> plt.show()  # doctest: +SKIP

.. image:: helioviewer-2.png

Next, let's see how we can download a JPEG 2000 image and load it into a SunPy
Map object.

The overall syntax is similar to the *download_png* request, expect instead of
specifying a single string to indicate which layers to use, here we
can specify the values as separate keyword arguments.::

   >>> from sunpy.net.helioviewer import HelioviewerClient
   >>> import matplotlib.pyplot as plt
   >>> from astropy.units import Quantity
   >>> from sunpy.map import Map
   >>> hv = HelioviewerClient()  # doctest: +REMOTE_DATA
   >>> data_sources = hv.get_data_sources()  # doctest: +REMOTE_DATA
   >>> filepath = hv.download_jp2('2012/07/05 00:30:00', observatory='SDO', instrument='HMI', detector=None, measurement='continuum')  # doctest: +REMOTE_DATA
   >>> hmi = Map(filepath)  # doctest: +REMOTE_DATA
   >>> xrange = Quantity([200, 550], 'arcsec')  # doctest: +REMOTE_DATA
   >>> yrange = Quantity([-400, 200], 'arcsec')  # doctest: +REMOTE_DATA
   >>> hmi.submap(xrange, yrange).peek()  # doctest: +SKIP

.. image:: helioviewer-3.png

Every JP2 file provided by the Helioviewer Project has been processed to generate an image that
can be used for browse purposes.  This typically involves following the standard image processing
procedure used by each instrument team to convert their science data into an image for a webpage.
The JP2 image is then scaled between 0 and 255 (byte-scaled).  Please note that the JP2 image data
is NOT the same as the original science data.  In the example above, SunPy queries Helioviewer for
the relevant JP2 file closest to the input time, downloads it, and selects a color table based on
the JP2 image meta data for plotting.  The color table is that used by the Helioviewer Project to
display JP2 images in their browse clients.

For more information about using querying Helioviewer.org, see the Helioviewer.org
API documentation at: `https://api.helioviewer.org/docs/v2/ <https://api.helioviewer.org/docs/v2/>`__.