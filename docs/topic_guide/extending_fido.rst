.. _sunpy-topic-guide-new-source-for-fido:

*******************************
Adding new data sources to Fido
*******************************

sunpy's data search and retrieval tool (``Fido``) is designed to be extensible, so that new sources of data or metadata can be supported, either inside or outside the sunpy core package.

There are two ways of defining a new client, depending on the complexity of the web service.
A "scraper" client inherits from `~sunpy.net.dataretriever.client.GenericClient` which provides helper methods for downloading from a list of URLs.
If the service you want to add has easily accessible HTTP or FTP URLs that have a well defined folder and filename structure, this is probably the best approach.
If the service you want to add requires making requests to an API with parameters for the search and getting a list of results in return, then you probably want to write a "full" client.

Before writing a new client, ensure you are familiar with how searches are specified by the `sunpy.net.attr` system, including combining them with logical operations.
When choosing a name for your new client it should have the form ``<name>Client`` as sunpy will split the name of the class to extract the name of your client.
The main place this is done is when constructing a `~.UnifiedResponse` object, where the name part can be used to index the response object.

.. _sunpy-topic-guide-new-source-for-fido-add-new-scraper-client:

Writing a new "scraper" client
==============================

A "scraper" Fido client (also sometimes referred to as a "data retriever" client) is a Fido client which uses the URL `~sunpy.net.scraper.Scraper` to find files on remote servers.
If the data provider you want to integrate does not provide a tree of files with predictable URLs then a "full" client is more likely to provide the functionality you need.

A new "scraper" client inherits from `~sunpy.net.dataretriever.client.GenericClient` and requires a minimum of these three components:

* A class method :meth:`~sunpy.net.base_client.BaseClient.register_values`; this registers the "attrs" that are supported by the client.
  It returns a dictionary where keys are the supported attrs and values are lists of tuples.
  Each `tuple` contains the "attr" value and its description.
* A class attribute ``baseurl``; this is a regular expression which is used to match the URLs supported by the client.
* A class attribute ``pattern``; this is a template used to extract the metadata from URLs matched by ``baseurl``.
  The extraction uses the `~sunpy.extern.parse.parse` format.

For a simple example of a scraper client, we can look at the implementation of `sunpy.net.dataretriever.sources.eve.EVEClient` in sunpy.
A version without documentation strings is reproduced below:

.. code-block:: python

    class EVEClient(GenericClient):
        baseurl = (r'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/'
                    r'L0CS/SpWx/%Y/%Y%m%d_EVE_L0CS_DIODES_1m.txt')
        pattern = '{}/SpWx/{:4d}/{year:4d}{month:2d}{day:2d}_EVE_L{Level:1d}{}'

        @classmethod
        def register_values(cls):
            from sunpy.net import attrs
            adict = {attrs.Instrument: [('EVE', 'Extreme ultraviolet Variability Experiment, which is part of the NASA Solar Dynamics Observatory mission.')],
                    attrs.Physobs: [('irradiance', 'the flux of radiant energy per unit area.')],
                    attrs.Source: [('SDO', 'The Solar Dynamics Observatory.')],
                    attrs.Provider: [('LASP', 'The Laboratory for Atmospheric and Space Physics.')],
                    attrs.Level: [('0', 'EVE: The specific EVE client can only return Level 0C data. Any other number will use the VSO Client.')]}
            return adict

This client scrapes all the URLs available under the base url ``http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/``.
`~sunpy.net.scraper.Scraper` is primarily focused on URL parsing based on time ranges, so the rest of the ``baseurl`` pattern specifies where in the pattern the time information is located, using `strptime <https://strftime.org/>`__ notation.
The ``pattern`` attribute is used to populate the results table from the URLs matched with the ``baseurl``.
It includes some of the time definitions, as well as names of attrs (in this case "Level").
The supported time keys are: 'year', 'month', 'day', 'hour', 'minute', 'second', 'millisecond'.

The attrs returned in the ``register_values()`` method are used to match your client to a search, as well as adding their values to the attr.
This means that after this client has been imported, running ``print(a.Provider)`` will show that the ``EVEClient`` has registered a provider value of ``LASP``.
In addition to this, a sanitized, lower cased version of the value will be available for tab completing, e.g. ``a.Provider.lasp`` or ``a.Level.zero``.

More Complex Clients
--------------------

Sometimes the attr values may not exist identically in the required URLs, and therefore can not be simply extracted with ``pattern``.
Say, for example, the Wavelength of a file is expressed in the URL as a passband by name; in this case conversion of the `~astropy.units.Quantity` object to the pass band name would be needed.
This is done addressed with the two following methods:

* :meth:`~sunpy.net.dataretriever.client.GenericClient.pre_search_hook` which will convert the passed attrs to their representation in the URL.
* :meth:`~sunpy.net.dataretriever.client.GenericClient.post_search_hook` which converts the retrieved metadata from a URL to the form in which they are desired to be represented in the response table.

A good example of the use of these two methods is the `sunpy.net.dataretriever.sources.norh.NoRHClient` in sunpy.

It may also be possible that the ``baseurl`` property needs to be customized based on attrs other than Time.
Since `~sunpy.net.scraper.Scraper` doesn't currently support generating directories that have non-time variables, the :meth:`~sunpy.net.dataretriever.client.GenericClient.search` needs to be customized.
The search method should in this case, generate a ``baseurl`` dependent on the values of these attrs, and then call ``super().search`` or `~sunpy.net.scraper.Scraper` for each ``baseurl`` generated.
For an example of a complex modification of the ``search()`` method see the implementation of `.SUVIClient.search`.

Customizing the Downloader
--------------------------

There is no method for a client creator to override the `parfive.Downloader` that is used to fetch the files.
This is because all downloads made by a single call to ``Fido.fetch`` share one instance of `parfive.Downloader`.
However, it is possible to pass keywords :meth:`parfive.Downloader.enqueue_file`, which is important if there is a need to customize the requests to a remote server, such as setting custom HTTP headers.
This is done by setting the ``enqueue_file_kwargs`` attribute of the client class.
One example from the `sunpy.net.dataretriever.sources.noaa.SRSClient` is:

.. code-block:: python

    class SRSClient(GenericClient):
        ...
        # Server does not support the normal aioftp passive command.
        enqueue_file_kwargs = {"passive_commands": ["pasv"]}
        ...

These keywords are passed to each call to :meth:`parfive.Downloader.enqueue_file`, so they will affect all files that are added for download by your client.

Examples
--------

Suppose any file of a data archive can be described by this URL ``https://some-domain.com/%Y/%m/%d/satname_{SatelliteNumber}_{Level}_%y%m%d%H%M%S_{any-2-digit-number}.fits``:

``baseurl`` becomes ``r'https://some-domain.com/%Y/%m/%d/satname_(\d){2}_(\d){1}_(\d){12}_(\d){2}\.fits'``.

Note all variables in the filename are converted to regex that will match any possible value for it.
A character enclosed within ``()`` followed by a number enclosed within ``{}`` is used to match the specified number of occurrences of that special sequence.
For example, ``%y%m%d%H%M%S`` is a twelve digit variable (with 2 digits for each item) and thus represented by ``r'(\d){12}'``.
Note that ``\`` is used to escape the special character ``.``.

``pattern`` becomes ``'{}/{year:4d}/{month:2d}{day:2d}/satname_{SatelliteNumber:2d}_{Level:1d}_{:6d}{hour:2d}{minute:2d}{second:2d}_{:2d}.fits'``.
Note the sole purpose of ``pattern`` is to extract the information from matched URL, using `~sunpy.extern.parse.parse`.
So the desired key names for returned dictionary should be written in the ``pattern`` within ``{}``, and they should match with the ``attr.__name__``.

``register_values()`` can be written as:

.. code-block:: python

    @classmethod
    def register_values(cls):

        from sunpy.net import attrs
        adict = {
        attrs.Instrument: [("SatName", "The description of Instrument")],
        attrs.Physobs: [('some_physobs', 'Phsyobs description')],
        attrs.Source: [('some_source', 'Source description')],
        attrs.Provider: [('some_provider', 'Provider description')],
        attrs.Level: [("1", "Level 1 data"), ("2", "Level 2 data")],
        attrs.SatelliteNumber: [("16", "Describe it"), ("17", "Describe it")]
        }

        return adict

.. _sunpy-topic-guide-new-source-for-fido-add-new-full-client:

Writing a "full" client
=======================

In this section we will describe how to build a "full" Fido client.
You should write a new "full" client if the data you are accessing can not be accessed via a URL template, for instance if you hit a web API with a query to return results for a search.

A new Fido client contains three major components:

* A subclass of `~sunpy.net.base_client.BaseClient` which implements ``search``, ``fetch``, and ``_can_handle_query``.
* Zero or more new `~sunpy.net.attr.Attr` classes to specify search parameters unique to your data source.
* An instance of `~sunpy.net.attr.AttrWalker` which can be used to walk the tree of `~sunpy.net.attr.Attr` instances and convert them into a form useful to your client's search method.

Search Attrs
------------

As described in `~sunpy.net.attr` the attr system allows the construction of complex queries by the user.
To make these complex queries easily processable by the clients the ``AttrWalker`` converts these into a set of queries which can be processed separately.
It does this by converting the input query to a set of queries which are ORed, but are complete queries.
This means the list of queries is an **OR** of **ANDs** (technically called `disjunctive normal form <https://en.wikipedia.org/wiki/Disjunctive_normal_form>`__).

Each query in the list of ORs contains all the information about that query so for example if the user provided a query like

.. code-block:: python

    a.Time("2020/02/02", "2020/02/03") & (a.Instrument("AIA") | a.Instrument("HMI"))

it would be passed to the client as

.. code-block:: python

    (a.Time("2020/02/02", "2020/02/03") & a.Instrument("HMI")) | (a.Time("2020/02/02", "2020/02/03") & a.Instrument("AIA"))

So you can process each element of the OR in turn without having to consult any other part of the query.

If the query the user provided contains an OR statement you get passed an instance of `~sunpy.net.attrs.AttrOr` and each sub-element of that `~sunpy.net.attrs.AttrOr` will be `~sunpy.net.attrs.AttrAnd` (or a single other attr class).
If the user query doesn't contain an OR you get a single `~sunpy.net.attr.Attr` instance or an `~sunpy.net.attrs.AttrAnd`.

For example you could get any of the following queries (using ``&`` for AND and ``|`` for OR):

* ``(a.Instrument("AIA") & a.Time("2020/02/02", "2020/02/03")) | (a.Instrument("HMI") & a.Time("2020/02/02", "2020/02/03"))``
* ``a.Time("2020/02/02", "2020/02/03")``
* ``a.Instrument("AIA") & a.Time("2020/02/02", "2020/02/03")``
* ``(a.Time(..) & a.Instrument("AIA") & a.Wavelength(30*u.nm, 31*u.nm)) | (a.Time(..) & a.Instrument("AIA") & a.Wavelength(30*u.nm, 31*u.nm))``

but you **would not** be passed queries which look like the following examples, even if that's how the user specified them:

* ``a.Time("2020/02/02", "2020/02/03") & (a.Instrument("AIA") | a.Instrument("HMI"))``
* ``a.Time(..) & (a.Instrument("AIA") | a.Instrument("AIA")) & a.Wavelength(30*u.nm, 31*u.nm))``

The Attr Walker
###############

Given the potential complexity of these combined attrs, converting them into other forms, such as query parameters or JSON etc involves walking the tree and converting each attr to the expected format in a given way.
This parsing and conversion of the query tree is deliberately not done using methods or attributes of the attrs themselves.
The attrs should be independent of any client in their implementation, so they can be shared between the different ``Fido`` clients.

A class is provided to facilitate this conversion, `~sunpy.net.attr.AttrWalker`.
The `~sunpy.net.attr.AttrWalker` class consists of three main components:

* **Creators**: The `~sunpy.net.attr.AttrWalker.create` method is one of two generic functions for which a different function is called for each Attr type.
  The intended use for creators is to return a new object dependent on different attrs.
  It is commonly used to dispatch on `~sunpy.net.attrs.AttrAnd` and `~sunpy.net.attrs.AttrOr`.

* **Appliers**: The `~sunpy.net.attr.AttrWalker.apply` method is the same as `~sunpy.net.attr.AttrWalker.create` in that it is a generic function.
  The only difference between it and `~sunpy.net.attr.AttrWalker.create` is its intended use.
  Appliers are generally used to modify an object returned by a creator with the values or information contained in other Attrs.

* **Converters**: Adding a converter to the walker adds the function to both the creator and the applier.
  For the VSO client this is used to convert each supported attr into a `~sunpy.net.attr.ValueAttr` which is then later processed by the appliers and creators.
  This pattern can be useful if you would otherwise have to repeat a lot of logic in each of the applier functions for each type of Attr you support.

An Example of ``AttrWalker``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example we will write a parser for some simple queries which uses `~sunpy.net.attr.AttrWalker` to convert the query to a `dict` of URL query parameters for a HTTP GET request.
Let's imagine we have a web service which you can do a HTTP GET request to ``https://sfsi.sunpy.org/search`` for some imaginary data from an instrument called SFSI (Sunpy Fake Solar Instrument).
This GET request takes three query parameters ``startTime``, ``endTime`` and ``level``, so a request might look something like: ``https://sfsi.sunpy.org/search?startTime=2020-01-02T00:00:00&endTime=2020-01-02T00:00:00&level=1``.
Which would search for level one data between 2020-01-01 and 2020-01-02.

As `~sunpy.net.attrs` has `~sunpy.net.attrs.Time` and `~sunpy.net.attrs.Level` we do not need to define any of our own attrs for this client.
We do however want to write our own walker to convert them to the form out client's ``search()`` method wants to send them to the server.

The first step is to setup the walker and define a creator method which will return a list of dicts, one for each independent search.

.. code-block:: python

    import sunpy.net.attrs as a
    from sunpy.net.attr import AttrWalker, AttrAnd, AttrOr, DataAttr

    walker = AttrWalker()

    @walker.add_creator(AttrOr)
    def create_or(wlk, tree):
        results = []
        for sub in tree.attrs:
            results.append(wlk.create(sub))

        return results

    @walker.add_creator(AttrAnd, DataAttr)
    def create_and(wlk, tree):
        result = dict()
        wlk.apply(tree, result)
        return [result]


The call ``wlk.apply(...)`` inside the creator will walk any nested attrs and add their values to the dictionary as defined by the applier registered to each attr type.
If we want our client to support searching by ``a.Time`` and ``a.Level`` as in the URL example above, we would need to register an applier for each of these attrs.

.. code-block:: python

    @walker.add_applier(a.Time)
    def _(wlk, attr, params):
        return params.update({'startTime': attr.start.isot,
                              'endTime': attr.end.isot})

    @walker.add_applier(a.Level)
    def _(wlk, attr, params):
        return params.update({'level': attr.value})


This combination of creators and appliers would allow support of any combination of queries consisting of ``a.Time`` and ``a.Level``.
Obviously, most clients would want to support more attrs than these two, and this could be done by adding more applier functions.

Adding "Attrs" to Registry
##########################

Registering of "attrs" ensures discoverability of search attributes supported by the corresponding sunpy Client.
For adding them to the Registry, we need to define a class method :meth:`~sunpy.net.base_client.BaseClient.register_values` that returns a dictionary of registered values.
This dictionary should have `~sunpy.net.attr.Attr` classes as keys and a list of tuples corresponding to that key representing the possible values the key "attr" can take.
Each tuple comprises of two elements.
The first one is a value and the second element contains a brief description of that value.
An example of writing ``register_values()`` for `~sunpy.net.dataretriever.client.GenericClient` is provided above.
Please note that it can be defined in a similar way for full clients too.

An Example of ``register_values()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    @classmethod
    def register_values(cls):

        from sunpy.net import attrs
        adict = {
        attrs.Instrument: [("LASCO", "Large Angle and Spectrometric Coronagraph")],
        attrs.Source: [('SOHO', 'Solar and Heliospheric Observatory')],
        attrs.Provider: [('SDAC', 'Solar Data Analysis Center')],
        attrs.Detector: [('C1', 'Coronograph 1'),
                         ('C2', 'Coronograph 2'),
                         ('C3', 'Coronograph 3')]
        }

        return adict

Registering custom attrs in the ``attrs`` namespace
---------------------------------------------------

When you have custom attrs defined in a separate attrs module, you can add them to the namespace using the :meth:`~sunpy.net.BaseClient._attrs_module` class method.
The method returns a tuple of length 2, where the first element is the target module name under which you want to add the custom attrs to the main attrs namespace.
The second is the import path to the source module where the custom attrs are defined.
Note that the source module here need not be an internal ``sunpy`` module, it could very well be external.
An example for this can be seen as implemented in the JSOC client:

.. code-block:: python

    @classmethod
    def _attrs_module(cls):
        return 'jsoc', 'sunpy.net.jsoc.attrs'

This adds all attrs that exist within ``sunpy.net.jsoc.attrs``, such as ``Keyword``, to ``attrs.jsoc``.
These can now be accessed via an import of the main attrs module, e.g., at ``a.jsoc.Keyword``.

Writing a Search Method
-----------------------

The ``search()`` method has the job of taking a set of user queries and returning an instance of `.QueryResponseTable` containing the results.

The general flow of a ``search()`` method is:

* Call your instance of an `.AttrWalker` to convert the input into a form expected by your API.
* Make as many requests to your API as needed to fulfill the query.
  Generally one per element of the outer `sunpy.net.attrs.AttrOr`.
* Process the response from your API into an instance of `.QueryResponseTable`.

To process the query with the `.AttrWalker`, call the :meth:`.AttrWalker.create` method:

.. code-block:: python

    def search(self, query):
        queries = walker.create(query)

Assuming the walker is the one we defined above, queries would be a list of dicts with the attrs processed into query parameters for the API URL.

.. note::

    If you want your search method to be able to be called independently of Fido, then you should accept a variable number of positional arguments (``*args``) and they should have the AND operator applied to them.
    This looks like:

    .. code-block:: python

        def search(self, *args):
            query = attr.and_(args)
            queries = walker.create(query)

Once the walker has processed the query into a form designed to be passed to your API, your ``search()`` method then needs to iterate over these parameters, make the requests, and process the results into a table.

In the following example we pretend our client has a method ``_make_search(query_parameters)`` which takes the query parameters and makes a request to our API.
We also pretend that the response is a json object in the form of a Python dictionary, which we want to put into the table.

.. code-block:: python

    def search(self, query):
        queries = walker.create(query)

        results = []
        for query_parameters in queries:
            results.append(self._make_search(query_parameters))

        return QueryResponseTable(results, client=self)

In reality, you probably want to post-process the results from your API before you put them in the table, they should be human readable first, with spaces and capitalization as appropriate.

Supporting file size estimates
##############################

The base client has a method for automatically estimating the total size of files in a given query: :meth:`~sunpy.net.base_client.QueryResponseTable.total_size`.
To enable to support for this, make sure the table returned by ``search`` has a column that contains filesizes as astropy quantities convertible to ``u.byte``, and set the ``size_column`` class attribute to the name of this column.

The ``_can_handle_query`` method
---------------------------------

The next required method is ``_can_handle_query``, this method tells ``Fido`` if your client might be able to return results for a given query.
If this method returns `True`, your clients ``search()`` method will be called for that query.
This method gets passed each query (in its independent form), and must either return ``True`` or ``False``.

A simple example, which just checks the type of ``attrs`` and not their values would be

.. code-block:: python

    @classmethod
    def _can_handle_query(cls, *query):
        query_attrs = set(type(x) for x in query)
        supported_attrs = {a.Time, a.Level}
        return supported_attrs.issuperset(query_attrs)

Note, that this method is a class method, it gets called without instantiating your client to speed up the dispatching.
If you are using the `~sunpy.net.dataretriever.client.GenericClient` as a base class, you do not need to implement this method, as it is already implemented in the base class.

Writing a Fetch Method
----------------------

The ``fetch()`` method of a Fido client is responsible for converting a set of search results (possibly sliced by the user) into a set of URLs to be downloaded.
Due to the history of clients and how they were implemented in sunpy, some existing clients support use outside of the ``Fido`` wrapper, this makes them appear more complex.
In this example we are going to write a ``fetch()`` method which is designed only to be called from ``Fido``.

The parameters for such a method should be:

.. code-block:: python

    def fetch(self, query_results, *, path, downloader, **kwargs):
    ...

The parameters here are:

* ``query_results`` which is an instance of `~.QueryResponseTable` or `~sunpy.net.base_client.QueryResponseRow`, these are the results the user wants to download.
* ``path=`` This is the path that the user wants the file to be downloaded to, this can be a template string (i.e. expects to have ``.format()`` called on it).
* ``downloader=`` This is a `parfive.Downloader` object which should be mutated by the ``fetch()`` method.
* ``**kwargs`` It is very important that ``fetch()`` methods accept extra keyword arguments that they don't use, as the user might be passing them to other clients via ``Fido``.

Processing the ``query_results`` Argument
#########################################

The ``query_results`` argument can be of two types `~.QueryResponseTable` or `~sunpy.net.base_client.QueryResponseRow`, as the user can slice the results table down to a single row and then pass that to ``Fido.fetch()``.
If you do not wish to handle a single row any differently to a table, you can place the `~sunpy.net.base_client.convert_row_to_table` decorator on your ``fetch()`` method which will convert the argument to a length one table when it is a single row object.

The primary function of the ``fetch()`` method is for you to convert this results object into a set of URLs for Fido to download.
This logic will be specific to your client.

Formatting the ``path=`` Argument
#################################

The path argument may contain format sections which are processed column names from the response table.
In addition to these it may contain the ``{file}`` format segment which is a placeholder for the filename.
Each row of the results table has a `~sunpy.net.base_client.QueryResponseRow.response_block_map` property which is a dictionary of valid format keys to values for that row.

In addition to the `~sunpy.net.base_client.QueryResponseRow.response_block_map` your fetch method also needs to be able to generate a filename for the file.
The simplest (but unlikely) scenario is that you know the filename for each file you are going to download before you do so, in this situation you would be able to generate the full filepath for each row of the response as follows

.. code-block:: python

    for row in query_results:
        filename = self._calculate_filename(row)
        filepath = path.format(file=filename, **row.response_block_map)

In the situation where you wish to be told the filename by the web server you are downloading the file from, it is a little more complex, you need to pass a callback function to :meth:`parfive.Downloader.enqueue_file` which will calculate the full filename in the context of the download, where the headers can be inspected for the filename the web server provides.

The filename callback passed to :meth:`parfive.Downloader.enqueue_file` accepts two arguments ``resp`` and ``url``.
``resp`` is an `aiohttp.ClientResponse` object which is returned when `parfive` requests the URL.
This response object allows us to inspect the headers of the response before the data is downloaded.
``url`` is the URL that was requested to generate the ``resp`` response.

To combine the formatting of the row with the extraction of the filename from the headers it is common to use `functools.partial` to generate many functions with different fixed parameters.
In the following example we will define a function which takes 4 arguments which we will use to generate the filename for the row.
This function will be called by `parfive` with the ``resp`` and ``url`` arguments.

.. code-block:: python

    def make_filename(path, row, resp, url):
        # Define a fallback filename based on the information in the search results
        name = f"row['ID'].fits"

        if resp:
            cdheader = resp.headers.get("Content-Disposition", None)
            if cdheader:
            _, params = sunpy.util.net.parse_header(cdheader)
            name = params.get('filename', "")

        return path.format(file=name, **row.response_block_map)

To reduce this function down to the two arguments expected we pre-specify the first two of these with `~functools.partial` before passing the function to `~parfive.Downloader.enqueue_file` inside the ``fetch()`` method.
Our simple example above now becomes:

.. code-block:: python

    for row in query_results:
        filepath = partial(make_filename, path, row)

Where the ``path`` variable is a `pathlib.Path` object provided as the ``path`` argument to ``fetch()``.

Adding URLs to be Downloaded
############################

For each file you wish for ``Fido`` to download (normally one per row of the ``query_results``) you need to call the :meth:`parfive.Downloader.enqueue_file` of the ``downloader`` argument.
Combining this with the simple example above it may look something like

.. code-block:: python

    for row in query_results:
        filename = self._calculate_filename(row)
        filepath = path.format(file=filename, **row.response_block_map)

        url = self._calculate_url(row)
        downloader.enqueue_file(url, filename=filepath)

If your filepath is a callback function, pass this to the ``filename=`` argument.

Your fetch method does not need to return anything, as long as ``enqueue_file`` is called for every file you want ``Fido`` to download.

Putting it all together
-----------------------

An example client class may look something like

.. code-block:: python

    import sunpy.util.net

    import sunpy.net.atrrs as a
    from sunpy.net.attr import AttrWalker, AttrAnd, AttrOr, DataAttr
    from sunpy.base_client import QueryResponseTable

    walker = AttrWalker()

    @walker.add_creator(AttrOr)
    def create_or(wlk, tree):
        results = []
        for sub in tree.attrs:
            results.append(wlk.create(sub))

        return results


    @walker.add_creator(AttrAnd, DataAttr)
    def create_and(wlk, tree):
        result = dict()
        wlk.apply(tree, result)
        return [result]


    @walker.add_applier(a.Time)
    def _(wlk, attr, params):
        return params.update({'startTime': attr.start.isot,
                                'endTime': attr.end.isot})


    @walker.add_applier(a.Level)
    def _(wlk, attr, params):
        return params.update({'level': attr.value})


    class ExampleClient(BaseClient):
        size_column = 'Filesize'

        def search(self, query):
            queries = walker.create(query)

            results = []
            for query_parameters in queries:
                results.append(self._make_search(query_parameters))

            return QueryResponseTable(results, client=self)

        def _make_filename(path, row, resp, url):
            # Define a fallback filename based on the information in the search results
            name = f"row['ID'].fits"

            if resp:
                cdheader = resp.headers.get("Content-Disposition", None)
                if cdheader:
                _, params = sunpy.util.net.parse_header(cdheader)
                name = params.get('filename', "")

            return path.format(file=name, **row.response_block_map)

        @convert_row_to_table
        def fetch(self, query_results, *, path, downloader, **kwargs):
            for row in query_results:
                filepath = partial(self._make_filename, path, row)

                url = f"https://sfsi.sunpy.org/download/{row['ID']}"
                downloader.enqueue_file(url, filename=filepath)

        @classmethod
        def _can_handle_query(cls, *query):
            query_attrs = set(type(x) for x in query)
            supported_attrs = {a.Time, a.Level}
            return supported_attrs.issuperset(query_attrs)
