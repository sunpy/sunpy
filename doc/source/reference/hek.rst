.. _hek:

===
HEK
===

Querying
^^^^^^^^
The HEK can be queried in a similar fashion as the VSO. Queries are
constructed from fundamental attributes using basic logic operations.
What the user ends up specifying is a logic condition the data must match
to be included in the result.

The primary interface to the HEK is the HEKClient, of which an object
must be created in order to use the HEK.

    >>> from sunpy.net import hek
    >>> client = hek.HEKClient()

Now we can use the client's query method in order to retrieve data from
the HEK. This example retrieves all AR events between the two specified
times. If query receives multiple arguments those are combined using AND.

    >>> res = client.query(hek.attrs.Time('2011-09-22T09:00:00', '2011-09-22T11:00:00'), hek.attrs.AR)

If we, for example, also wanted CMEs in this period.


    >>> res = client.query(hek.attrs.Time('2011-09-22T09:00:00', '2011-09-22T11:00:00'), hek.attrs.AR | hek.attrs.CE)

Which reads out as "find datasets that are within the given dataset AND
are (an AR OR a CE)".

Module documentation
^^^^^^^^^^^^^^^^^^^^

.. currentmodule:: sunpy.net.hek

.. automodule:: sunpy.net.hek
   :members:
