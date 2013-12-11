------------------------
Using SunPy's VSO module
------------------------

The `Virtual Solar Observatory (VSO) <http://virtualsolar.org>`_ 
is a service which presents a
homogenoeous interface to heterogeneous data-sets and services.  Using
the VSO, a user can query multiple data providers simultaneously, and
then download the relevant data.  SunPy uses the VSO through the ``vso``
module, which was developed through support from the `European Space
Agency Summer of Code in Space (ESA-SOCIS) 2011 
<http://sophia.estec.esa.int/socis2011/>`_.

Setting up the VSO interaction
---------------------------------

SunPy's VSO module is in ``sunpy.net``.  It can be imported into your
Python session as follows:

    >>> from sunpy.net import vso
    >>> client=vso.VSOClient()

Obtaining data via the VSO is essentially a two-stage process.  In the
first stage, you ask the VSO to find the data you want.  The VSO
queries various data-providers looking for your data.  In the second
stage, you download the data, if there is any data that matches your
request.  The VSO client handles the particulars of how the data from
the data provider is downloaded to your computer.

Querying VSO
------------

VSO module can be used in two different ways, a new and powerful 
query style (sunpy's way) and a simpler one based on the 
Solarsoft/IDL's VSO query client, VSO_SEARCH.pro (legacy).

To search the VSO, you need at least start time, an end time, 
and an instrument. Different examples are shown in the following sections.

SunPy's way
^^^^^^^^^^^

Let's start with a very simple query.  We could ask for all SOHO/EIT 
data between January 1st and 2nd, 2001.

    >>> qr = client.query(vso.attrs.Time('2001/1/1', '2001/1/2'), vso.attrs.Instrument('eit'))

The variable ``qr`` is a Python list of
response objects, each one of which is a record found by the VSO. How
many records have been found?  You can find that out be typing

    >>> qr.num_records()
    122

To get a little bit more information about the records found, try

    >>> qr.show() # doctest:+SKIP
    ...


Now, let's break down the arguments of ``client.query`` to understand
better what we've done.  The first argument:

    ``vso.attrs.Time('2001/1/1', '2001/1/2')``

sets the start and end times for the query (any date/time 
format understood by SunPy's :ref:`parse_time function <parse-time>` 
can be used to specify dates and time).  The second argument:

    ``vso.attrs.Instrument('eit')``

sets the instrument we are looking for.  So what is going on here?
The notion is that a VSO query has a set of attribute objects -
described in ``vso.attrs`` - that are specifed to construct the query.
For the full list of vso attributes, use

    >>> help(vso.attrs) # doctest:+SKIP

Note that due to quirks at the VSO, we do not recommend that the
extent object ``vso.attrs.Extent`` be in your query.  Instead, we
recommend that any extent filtering you need to do be done on the
queries made without setting a value to the ``vso.attrs.Extent`` object.
As we will see, this query style can take more than two arguments,
each argument separated from the other by a comma.  Each of those
arguments are chained together using a logical ``AND``.

This query style allows you to combine these VSO attribute objects
in complex ways that are not possible with the legacy query style.

So, let's look for the EIT and MDI data on the same day:

    >>> qr=client.query(vso.attrs.Time('2001/1/1', '2001/1/2'), vso.attrs.Instrument('eit') | vso.attrs.Instrument('mdi'))
    >>> qr.num_records()
    144
    >>> qr.show() # doctest:+SKIP
    ...

The two instrument types are joined together by the operator "|".
This is the ``or`` operator.  Think of the above query as setting a set
of conditions which get passed to the VSO.  Let's say you want all the
EIT data from two separate days:

    >>> qr=client.query(vso.attrs.Time('2001/1/1', '2001/1/2') | vso.attrs.Time('2007/8/9', '2007/8/10'), vso.attrs.Instrument('eit') )
    >>> qr.num_records()
    227

Each of the arguments in this query style can be thought of as
setting conditions that the returned records must satisfy.  You can
set the wavelength; for example, to return the 171 Angstrom EIT results

    >>> qr=client.query(vso.attrs.Time('2001/1/1', '2001/1/2'), vso.attrs.Instrument('eit'), vso.attrs.Wave(171,171) )
    >>> qr.num_records()
    4

Legacy query (ssw)
^^^^^^^^^^^^^^^^^^

Let's show now how the first example is executed with the legacy mode.

So, as before,  we want EIT data between 2001/01/01 and 2001/01/02

    >>> qr=client.query_legacy(tstart='2001/01/01', tend='2001/01/02', instrument='EIT')

which is almost identical to what you would type in a Solarsoft/IDL
session.  So, what's happening with this command?  The client is going
out to the web to query the VSO to ask how many files EIT images are
in the archive between the start of 2001/01/01 and the start of
2001/01/02.  The same query can also be performed using a slightly different
syntax.  For example

    >>> qr=client.query_legacy('2001/1/1', '2001/1/2', instrument='EIT')

both gives the same result. The variable ``qr`` is a Python list of
response objects, each one of which is a record found by the VSO. How
many records have been found?  You can find that out be typing

    >>> qr.num_records()
    122

To get a little bit more information, try

    >>> qr.show() # doctest:+SKIP
    ...

The Solarsoft legacy query has more keywords available: to find out
more about the legacy query, type: 

    >>> help(client.query_legacy) # doctest:+SKIP

As an example, let's say you just want the EIT 171 Angstrom files for
that data.  These files can be found by

    >>> qr=client.query_legacy(tstart='2001/01/01', tend='2001/01/02', instrument='EIT', min_wave='171', max_wave='171', unit_wave='Angstrom')

which yields four results, the same as the VSO IDL client.

Downloading data
----------------
Both query styles return the same type of response.  
This means you can use the same command and syntax to download your data.

Having located the data you want, you can download it using the
following command:

    >>> res=client.get(qr, path='/ThisIs/MyPath/to/Data/{file}.fits')

This downloads the query results into the directory
``/ThisIs/MyPath/to/Data`` naming each downloaded file with the
filename ``{file}`` obtained from the VSO , and appended with the suffix
``.fits``.  The ``{file}`` option uses the file name obtained by the VSO
for each file.  You can also use other properties of the query return
to define the path where the data is saved.  For example, to save the
data to a subdirectory named after the instrument, use

    >>> res=client.get(qr, path='/ThisIs/MyPath/to/Data/{instrument}/{file}.fits')

Note that the download process is spawned in parallel to your existing
Python session.  This means that the remainder of your Python script
will continue as the download proceeds.  This may cause a problem if
the remainder of your script relies on the presence of the downloaded
data.  If you want to resume your script after all the data has been
downloaded then append ``.wait()`` to the ``get`` command above, i.e.,

     >>> res=client.get(qr, path='/Users/ireland/Desktop/Data/{instrument}/{file}.fits').wait()

More information on the options available can be found through the
standard Python ``help`` command.



Finally, please let us know if you encounter a bug while using the VSO
capabilities of SunPy.  Bugs are best reported via the issue tracking
system on GitHub - go to https://github.com/sunpy/sunpy/issues and
click on the "New Issue" button.
