The SunPy project is happy to announce the release of SunPy 0.6.0.
This is a major SunPy release with lots of changes that will make SunPy even
better than it was before.
This release consists of 1,139 commits from 24 different people and
9 new contributors!

The major changes in this release are:

    * Most functions throughout the SunPy code base expect Astropy
      Quantity objects, and return Astropy Quantity objects.
    * Python 2.6 support has ended, we do not expect this release to
      work under Python 2.6.
    * Sample data has been removed from SunPy but helpers for
      downloading sample data have been added to sunpy.data.
    * TimeRange has a new property based API, e.g. start and end are
      now properties.
    * Map.rotate() now pads arrays to prevent loss of data under
      rotation.
    * Map.rotate() now defaults to the slower but more accurate
      bi-quartic interpolation method (order=4).
    * SunPy colormaps are now registered with matplotlib, allowing
      their use from imshow and similar functions after the import
      of sunpy.cm.
    * JSOC client export calls now checks for an email address.
    * Solar rotation calculation functionality has been added, along
      with functionality to de-rotate MapCubes.
    * Files downloaded through the VSO now keep their file
      extension.

The people who have contributed to this release are:

    Stuart Mumford
    Daniel Ryan
    Jack Ireland
    Steven Christe
    Albert Y. Shih
    Asish Panda
    Andrew Inglis
    Rishabh Sharma
    David Perez-Suarez
    Rajul Srivastava
    Ruben De Visscher *
    Dumindu Buddhika *
    Andrew Leonard
    Goran Cetusic *
    Ishtyaq Habib *
    Nabil Freij
    Simon Liedtke
    Abigail Stevens *
    Ambar Mehrotra *
    Erik M. Bray *
    Jaylen Wimbish *
    Larry Manley
    Norbert Gyenge
    Rishabh Mishra *

Where a * indicates their first contribution.
