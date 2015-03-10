"""
This module is broken into two layers. At the bottom layer individual clients
(LightCurve clients, VSO) operate and download files.
At the top level the Factory instance decides which client/s can best serve
the query. The factory instance breaks up a query into modular pieces which could be served by
clients at the lower layer.  It then calls an instance of the required client
to download the data.  New individual clients can be registered with factory instance.
The clients should have _can_handle_query class method which returns a boolean
representing whether the client can handle modularised query.

Examples
--------
>>> from sunpy.net.unifieddownloader import UnifiedDownloader
>>> import sunpy.net.vso.attrs as attrs
>>> results = UnifiedDownloader.query(attrs.Time("2012/1/1", "2012/1/2"), attrs.Instrument('lyra'))

query method returns UnifiedResponse object. This is a subclass of List.
__str__() method has been overloaded to show all the files downloaded by multiple
clients. The numfile property shows the number of files that need to be downloaded.
Each element in the container is a QueryResponse object returned by the underneath
client. This QueryResponse object has client member holding an instance of client.

>>> print results
Start time  End time    Source  Instrument  URL
----------  --------    ------  ----------  ---
2012/01/01  2012/01/02  Proba2  lyra        http://proba2.oma.be/lyra/data/bsd/2012/01/01/lyra_20120101-000000_lev2_std.fits
2012/01/02  2012/01/03  Proba2  lyra        http://proba2.oma.be/lyra/data/bsd/2012/01/02/lyra_20120102-000000_lev2_std.fits
>>> print results.numfile
2
>>> print len(results)
1 #Indicating only a single client was used to service the query.

>>> downresp = Downloader.get(results)
>>> downresp.wait()

get method returns DownloadResponse object.This is list of Results object (same ones as
the VSO Results object). It has a wait method which returns a list of file paths after
completion of downloading of all files.

Notes
-----

1)In case multiple clients can serve the query, factory instance delegates the query to one of them.
Eg. EVE data can be sourced from VSO and EVEClient. Here EVEClient is delegated the duty as it is faster
at downloading data.

2) In case a Time Range is specified, files downloaded might span complete days. A Truncate method
from LightCurve factory should be used to get specific data.
"""

from ..attrs import *

from .downloader_factory import Fido

import clients