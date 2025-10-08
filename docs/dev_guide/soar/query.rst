.. _sunpy-soar-dev-guide-query:

***************************************
Request methods ``sunpy-soar`` supports
***************************************

sunpy-soar currently supports two REQUEST methods: ``doQuery`` and ``doQueryFilteredByDistance``.

``doQuery``: This is the standard method used when no specific distance attribute is included in the search query.
It performs a general query based on the provided parameters, retrieving the data that matches the criteria specified.

``doQueryFilteredByDistance``: This method is employed when a distance parameter is included in the search query.
Unlike ``doQuery``, this method filters the entire database based on the specified distance value.
The time attribute is not necessarily required when using ``doQueryFilteredByDistance``.
The distance range of values is appended to the end of the query using ``&DISTANCE(dmin, dmax)``, where ``dmin`` and ``dmax`` are Astropy quantities representing distances.
These values must fall within the range of 0.28 AU to 1.02 AU; otherwise, the query will not return any results.

Using the example below,

.. code-block:: python

    >>> import astropy.units as u
    >>> import sunpy.net.attrs as a
    >>> from sunpy.net import Fido

    >>> instrument = a.Instrument("RPW")
    >>> level = a.Level(2)
    >>> distance = a.soar.Distance(0.28 * u.AU, 0.30 * u.AU)
    >>> result = Fido.search(instrument & level & distance) # doctest: +REMOTE_DATA
    >>> result # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    680 Results from the SOARClient:
    <BLANKLINE>
    Instrument     Data product    Level        Start time               End time        Filesize SOOP Name
                                                                                          Mbyte
    ---------- ------------------- ----- ----------------------- ----------------------- -------- ---------
           RPW rpw-tds-surv-hist2d    L2 2022-10-09 00:00:00.000 2022-10-10 00:00:00.000    0.084      none
    ...
    <BLANKLINE>
    <BLANKLINE>

Here the query's "REQUEST" type to "doQueryFilteredByDistance", which is a special method that filters the entire database based on the specified distance value.

The actual query this example produces is,

.. code-block:: python

   "SELECT+*+FROM+v_sc_data_item+WHERE+instrument='RPW'+AND+level='L2'&DISTANCE(0.28,0.30)"

How can other request methods be added?
=======================================

To add support for additional request methods, you'll need to consider the impact they will have on the query structure.
The REQUEST parameter in the query must be updated accordingly, and any necessary modifications to the query string should be made to accommodate the new method.
If the new request method requires a specific attribute, you may need to add it as a class in the attrs.py file (if it doesn't already exist in (if it doesn't already exist in `sunpy.net.attrs <https://github.com/sunpy/sunpy/blob/main/sunpy/net/attrs.py/>`__).
Additionally, a walker for this attribute will need to be created.
