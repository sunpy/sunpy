---------------------------------------
Finding and Downloading Data using Fido
---------------------------------------

This guide outlines how to search for and download data using SunPy's
Federated Internet Data Obtainer...or more usually (and
sanely) referred to as Fido.  Fido is a unified interface for seeking
and downloading solar physics data irrespective of the underlining
client or webservice through which the data is obtained, e.g. VSO,
JSOC, etc.  It therefore supplies a single, easy and consistent way to
obtain most forms of solar physics data.

Import
------

SunPy's Fido module is in ``sunpy.net``.  It can be imported as follows::

    >>> from sunpy.net import Fido, attrs as a

Searching for Data Using Fido
-----------------------------

To search for data with Fido, your query needs at minimum a start time,
an end time, and an instrument.  Enter these properties using SunPy's
attrs module::

    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument('lyra'))

This returns an `sunpy.net.fido_factory.UnifiedResponse` object
containing information on the available online files which fit the
criteria specified by the attrs objects in the above call.  It does
not download the files.  For instructions on how to download data
using Fido, see :ref:`downloading_data`. 

To see a summary of results of our query, simple type the name of the
variable set to the Fido search, in this case, result::

    >>> result
    
Queries can be made more flexible or specific by adding more attrs
objects to the Fido search.  Specific passbands taken with given
instruments can be searched for by supplying an astropy Quantity to
the Wavelength attribute::

    >>> import astropy.units as u
    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument('norh'), a.Wavelength(17*u.GHz))

Data of a given cadence can also be specified using the Sample
attribute::

    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument('aia'), a.Wavelength(171*u.angstrom), a.Sample(10*u.minute))

To search for data from multiple instruments, wavelengths, times etc.,
use the pipe ``|`` operator.  This joins queries together just as the
logical ``OR`` operator would::

    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument('lyra') | a.Instrument('rhessi'))

    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument('aia'), a.Wavelength(171*u.angstrom) | a.Wavelength(94*u.angstrom))

The above query searches the Virtual Solar Observatory (VSO) showing 
that Fido can be used instead of explicitly searching the VSO.

.. _downloading_data:

Downloading data
----------------
Once you have located your files via a ```Fido.search```, you can download
them via ```Fido.fetch```::

    >>> downloaded_files = Fido.fetch(result)

This downloads the files to the location set in you sunpy config
file.  It also returns a list ``downloaded_files``, of absolute file paths
of where the files have been downloaded to.
