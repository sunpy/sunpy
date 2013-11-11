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
the HEK.   The default location for the HEK is
*http://www.lmsal.com/hek/her*. This can be changed when the client is
created, for example,

    >>> client2 = hek.HEKClient(url='http://www.the_other_hek.org')

The following example retrieves all AR events between the two specified
times. If query receives multiple arguments those are combined using AND.

    >>> res = client.query(hek.attrs.Time('2011-09-22T09:00:00', '2011-09-22T11:00:00'), hek.attrs.AR)

If we, for example, also wanted CMEs in this period.

    >>> res = client.query(hek.attrs.Time('2011-09-22T09:00:00', '2011-09-22T11:00:00'), hek.attrs.AR | hek.attrs.CE)

Which reads out as "find datasets that are within the given dataset AND
are (an AR OR a CE)".

The result of an HEK query is a list of dictionaries.  Each dictionary
contains the information on that feature/event as held in the HEK.


Module documentation
^^^^^^^^^^^^^^^^^^^^

:class:'HEKClient'
=============

The main class that interacts with the HEK service.  Provides a method
that queries the HEK service, a download method, and a merge method
(this method ensures that a unique list of events is returned).


:class:'Response'
============

Handles the response from the HEK service.  It retrieves the VOEvent
object associated with a given event and returns it as either a Python
dictionary or an XML string.  It also defines some very simple
properties of each HEK feature/event: the start and end time of the
event, and the instrument that observed it.


HEK attributes
^^^^^^^^^^^^


.. autoclass: sunpy.net.hek.HEKClient

.. currentmodule:: sunpy.net.hek

.. automodule:: sunpy.net.hek
   :members:
