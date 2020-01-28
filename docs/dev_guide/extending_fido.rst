.. _extending_fido:

***************************************
Extending Fido with New Sources of Data
***************************************

The `~sunpy.net.fido_factory.UnifiedDownloaderFactory` (``Fido``) object is extensible with new clients, which can interface with webservices and download new data.
There are two ways of defining a new client, depending on the complexity of your webservice.
A "simple" client inherits from `sunpy.net.dataretriver.GenericClient` which provides helper methods for downloading from a list of URLs.
If your webservice provides a list of HTTP or FTP urls that can easily be obtained from a search, this is probably the route to go.
If your webservice requires you to do complex parsing of the search, or needs to construct specific objects to interface with the webservice, or you need control over the download implementation (i.e. does not just return a list of URLs) then you probably want to write a "full" client.


Writing a new "simple" client
=============================



Writing a full client
=====================


A new Fido client contains three major components:

* A subclass of `sunpy.net.base_client.BaseClient` which implements the interface defined on that `~abc.ABC`, namely ``search``, ``fetch``, and ``_is_datasource_for``.
* Zero or more new `sunpy.net.attr.Attr` classes to specify search parameters unique to your data source.
* An instance of `sunpy.net.attr.AttrWalker` which can be used to walk the tree of `~sunpy.net.attr.Attr` instances and convert them into a form useful to your client's search method.


Processing Search Attrs
-----------------------

As described in `sunpy.net.attr` the attr system allows the construction of complex queries by the user.
It then converts them to `disjuntive normal form <https://en.wikipedia.org/wiki/Disjunctive_normal_form>`__ an **OR** of **ANDS**.
This means that as a client author, when you get passed a query (which contains an OR statement), the outer most `~sunpy.net.attr.Attr` is `~sunpy.net.attr.AttrOr` and each sub-tree of the `~sunpy.net.attr.AttrOr` will be `~sunpy.net.attr.AttrAnd` (or a single other attr class).
For example you could get any of the following queries (using ``&`` for AND and ``|`` for OR):

* ``(a.Instrument("AIA") & a.Time("2020/02/02", "2020/02/03")) | (a.Instrument("HMI") & a.Time("2020/02/02", "2020/02/03"))``
* ``a.Time("2020/02/02", "2020/02/03")``
* ``a.Instrument("AIA") & a.Time("2020/02/02", "2020/02/03")``
* ``(a.Time(..) & a.Instrument("AIA") & a.Wavelength(30*u.nm, 31*u.nm)) | (a.Time(..) & a.Instrument("AIA") & a.Wavelength(30*u.nm, 31*u.nm))``

but you **would not** be passed queries which look like the following examples, even if that's how the user specified them:

* ``a.Time("2020/02/02", "2020/02/03") & (a.Instrument("AIA") | a.Instrument("HMI"))``
* ``a.Time(..) & (a.Instrument("AIA") | a.Instrument("AIA")) & a.Wavelength(30*u.nm, 31*u.nm))``

Presenting the query in this form makes it suitable for many different kinds of webservices, as each block which is ORed normally maps to a single search request.


The Attr Walker
###############

Given the potential complexity of these combined attrs, converting them into other forms, such as query parameters or JSON etc involves walking the tree and converting each attr to the expected format in a given way.
This parsing and conversion of the query tree is deliberately not done using methods or attributes of the attrs themselves.
The attrs should be independent of any client in their implementation, so they can be shared between the different ``Fido`` clients.

A class is provided to facilitate this conversion, `sunpy.net.attr.AttrWalker`.
The `~sunpy.net.attr.AttrWalker` class consists of three main components:

* **Creators**: The `~sunpy.net.attr.AttrWalker.create` method is one of two generic functions for which a different function is called for each Attr type.
  The intended use for creators is to return a new object dependant on different attrs.
  It is commonly used to dispatch on `~sunpy.net.attr.AttrAnd` and `~sunpy.net.attr.AttrOr`.

* **Appliers**: The `~sunpy.net.attr.AttrWalker.apply` method is the same as `~sunpy.net.attr.AttrWalker.create` in that it is a generic function.
  The only difference between it and `~sunpy.net.attr.AttrWalker.create` is its intended use.
  Appliers are generally used to modify an object returned by a creator with the values or information contained in other Attrs.

* **Converters**: Adding a converter to the walker, adds the function to the both the creator and the applier.
  For the VSO client this is used to convert each supported attr into a `sunpy.net.attr.ValueAttr` which is then later processed by the appliers and creators.
  This pattern can be useful if you would otherwise have to repeat a lot of logic in each of the applier functions for each type of Attr you support.


An Example of ``AttrWalker``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example we will write a parser for some simple queries which uses `~sunpy.net.attr.AttrWalker` to convert the query to a `dict` of URL query parameters for a HTTP GET request.
Let's imagine we have a webservice which you can do a HTTP GET request to ``https://sfsi.sunpy.org/search`` for some imaginary data from an instrument called SFSI (SunPy Fake Solar Instrument).
This GET request takes three query parameters ``startTime``, ``endTime`` and ``level``, so a request might look something like: ``https://sfsi.sunpy.org/search?startTime=2020-01-02T00:00:00&endTime=2020-01-02T00:00:00&level=1``.
Which would search for level one data between 2020-01-01 and 2020-01-02.

As `sunpy.net.attrs` has `~sunpy.net.attrs.Time` and `~sunpy.net.attrs.Level` we don't need to define any of our own attrs for this client.
We do however want to write our own walker to convert them to the form out client's ``search()`` method wants to send them to the server.


The first step is to setup the walker and define a creator method which will return either a dict (for a single query) or a list of dicts for multiple queries.

..code-block:: python

    import sunpy.net.atrrs as a
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
        return result
