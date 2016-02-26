"""
This module is broken into two layers. At the bottom layer individual clients
(LightCurve clients, VSO) operate and download files.
At the top level the Factory instance decides which client/s can best serve
the query. The factory instance breaks up a query into modular pieces which
could be served by clients at the lower layer.  It then calls an instance of
the required client to download the data.  New individual clients can be
registered with factory instance.
Each client should have _can_handle_query class method which returns a boolean
representing whether the client can handle a modularised query.

Examples
--------
>>> from sunpy.net import Fido, attrs as a
>>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"), a.Instrument('lyra'))

The query method returns UnifiedResponse object, which is a subclass of List.
The file_num property shows the number of files that need to be downloaded.
Each element in the container is a QueryResponse object returned by the
underneath client. This QueryResponse object has client member holding an
instance of client.

>>> print results
[<Table masked=False length=2>
     Start Time           End Time       Source  Instrument
     string152           string152      string48  string32
------------------- ------------------- -------- ----------
2012-01-01 00:00:00 2012-01-02 00:00:00   Proba2       lyra
2012-01-02 00:00:00 2012-01-03 00:00:00   Proba2       lyra]

>>> print results.file_num
2
>>> print len(results)
1 # Indicating only a single client was used to service the query.

>>> downresp = Fido.get(results)
>>> files = downresp.wait()

The get method returns a DownloadResponse object.  This is list of Results
object (same ones as the VSO Results object). It has a wait method which
returns a list of file paths after completion of downloading of all files.

Notes
-----

1) When multiple clients can serve the query, the factory instance delegates
the query to one of them. In the case when there are exactly two
clients that can serve the query and one of them is the VSO, and the other one
is a web service accessible via a SunPy instrument-specific client, the SunPy
instrument-specific client is delegated the duty to get the data.  This is
because it is generally faster to get data directly to instrument's web
service, rather than going indirectly through the VSO.

Consider the following example. EVE data is available either through the VSO
or by accessing directly the website that hold EVE data.  The SunPy EVEClient
connects to that website.  Therefore in this case, SunPy's EVEClient will be
used to service the query and get the data.

2) In case a Time Range is specified, files downloaded might span complete
days. A Truncate method from LightCurve factory should be used to get specific
data.
"""

from .client import QueryResponseBlock, QueryResponse, GenericClient
from .downloader_factory import DownloadResponse, UnifiedDownloaderFactory, Fido

from . import clients
