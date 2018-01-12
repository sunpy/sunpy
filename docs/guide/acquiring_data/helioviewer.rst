-----------------------------------
Querying Helioviewer.org with SunPy
-----------------------------------
SunPy can be used to make several basic requests using the The `Helioviewer.org API <https://helioviewer.org/api/>`_
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
the API using the same parameters one would use when making a web request.

Nearly all requests require the user to specify the data they are interested in
and this can be done using one of two methods:

1. Call "get_data_sources()" to get a list of the data that is available, and use the source id numbers referenced in the result to refer to a particular dataset, or,
2. Specify the four components of a Helioviewer.org data source or layer: *observatory*, *instrument*, *detector* and *measurement*.

Let's begin by getting a list of data sources available on the server
using the get_datasources method::

    >>> from sunpy.net.helioviewer import HelioviewerClient

    >>> hv = HelioviewerClient()
    >>> datasources = hv.get_data_sources()  # doctest: +REMOTE_DATA

    >>> # print a list of datasources and their associated ids
    >>> for observatory, instruments in datasources.items():  # doctest: +REMOTE_DATA
    ...     for inst, detectors in instruments.items():  # doctest: +REMOTE_DATA
    ...         for det, measurements in detectors.items():  # doctest: +REMOTE_DATA
    ...             for meas, params in measurements.items():  # doctest: +REMOTE_DATA
    ...                 print("%s %s: %d" % (observatory, params['nickname'], params['sourceId']))  # doctest: +REMOTE_DATA
    SOHO EIT 171: 0
    SOHO EIT 195: 1
    SOHO EIT 284: 2
    SOHO EIT 304: 3
    SOHO LASCO C2: 4
    SOHO LASCO C3: 5
    SOHO MDI Mag: 6
    SOHO MDI Int: 7
    SDO AIA 94: 8
    SDO AIA 131: 9
    SDO AIA 171: 10
    SDO AIA 193: 11
    SDO AIA 211: 12
    SDO AIA 304: 13
    SDO AIA 335: 14
    SDO AIA 1600: 15
    SDO AIA 1700: 16
    SDO AIA 4500: 17
    SDO HMI Int: 18
    SDO HMI Mag: 19
    STEREO_A EUVI-A 171: 20
    STEREO_A EUVI-A 195: 21
    STEREO_A EUVI-A 284: 22
    STEREO_A EUVI-A 304: 23
    STEREO_A COR1-A: 28
    STEREO_A COR2-A: 29
    STEREO_B EUVI-B 171: 24
    STEREO_B EUVI-B 195: 25
    STEREO_B EUVI-B 284: 26
    STEREO_B EUVI-B 304: 27
    STEREO_B COR1-B: 30
    STEREO_B COR2-B: 31
    PROBA2 SWAP 174: 32
    Yohkoh SXT AlMgMn: 33
    Yohkoh SXT thin-Al: 34
    Yohkoh SXT white-light: 35

At time of writing (2014/01/06) Helioviewer provides JP2 images from AIA, HMI, LASCO C2/C3, EIT,
MDI, STEREO A/B COR1/2 & EUVI, SWAP and SXT.  New sources of JP2 images are being added every few months;
please use the code snippet above to get an up-to-date list of available data sources.


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
   >>> filepath = hv.download_jp2('2012/07/05 00:30:00', observatory='SDO', instrument='HMI', detector='HMI', measurement='continuum')  # doctest: +REMOTE_DATA
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
API documentation at: `https://helioviewer.org/api/ <https://helioviewer.org/api/>`__.
