------------------------
Using SunPy's VSO module
------------------------

The Virtual Solar Observatory (VSO) is a service which presents a
homogenoeous interface to heterogeneous data-sets and services.  Using
the VSO, a user can query multiple data providers simultaneously, and
then download the relevant data.  SunPy uses the VSO through the 'vso'
module, which was developed through support from the European Space
Agency Summer of Code in Space (ESA-SOCIS) 2011.

1. Setting up the VSO interaction
--------------------------

SunPy's VSO module is in sunpy.net.  It can be imported into your
IPython session as follows:

    >>> from sunpy.net import vso


2. A simple query - using the legacy query
---------------------------------

Obtaining data via the VSO is essentially a two-stage process.  In the
first stage, you ask the VSO to find the data you want.  The VSO
queries various data-providers looking for your data.  In the second
stage, you download the data, if there is any data that matches your
request.  The VSO client handles the particulars of how the data from
the data provider is downloaded to your computer.

Let's start with a very simple query.  To search the VSO, you need a
start time, an end time, and an instrument. We have provided two
different syntaxes for doing this search.  The first query syntax
basically copies what you already may be used to from Solarsoft/IDL's
VSO query client, VSO_SEARCH.pro.  This is known as a 'legacy' query,
purely because the syntax is based on the legacy of Solarsoft/IDL's
VSO query client.  The second query syntax is much more powerful, but
slightly less familiar to Solarsoft/IDL users (which is why we have
two different syntaxes).

Let's say I want all the EIT data between 2001/01/01 and 2001/01/02
(note that you can use any format understood by the parse_time
function to specify the date and times you want).  Using the legacy
query syntax, this is simply

    >>> client=vso.VSOClient()
    >>> qr=client.query_legacy(tstart='2001/01/01', tend='2001/01/02', instrument='EIT')

which is almost identical to what you would type in a Solarsoft/IDL
session.  So, what's happening with this command?  The client is going
out to the web to query the VSO to ask how many files EIT images are
in the archive between the start of 2001/01/01 and the start of
2001/01/02.  The same query can also be performed using a slightly different
syntax.  For example

    >>> qr=client.query_legacy('2001/1/1', '2001/1/2', instrument='EIT')

both gives the same result. The variable "qr" is a Python list of
response objects, each one of which is a record found by the VSO. How
many records have been found?  You can find that out be typing

    >>> qr.num_records()
    122

To get a little bit more information, try

    >>> qr.show()

The Solarsoft legacy query has more keywords available: to find out
more about the legacy query, type: 

    >>> help(client.query_legacy)

As an example, let's say you just want the EIT 171 Angstrom files for
that data.  These files can be found by

    >>> qr=client.query_legacy(tstart='2001/01/01', tend='2001/01/02', instrument='EIT', min_wave='171', max_wave='171', unit_wave='Angstrom')

which yields four results, the same as the VSO IDL client.

3. Downloading the data
--------------------

Having located the data you want, you can download it using the
following command:

    >>> res=client.get(qr, path='/Users/ireland/Desktop/Data/{file}.fits')

This downloads the query results into the directory
/Users/ireland/Dekstop/Data naming each downloaded file with the
filename '{file}' obtained from the VSO , and appended with the suffix
'.fits'.  The '{file}' option uses the file name obtained by the VSO
for each file.  You can also use other properties of the query return
to define the path where the data is saved.  For example, to save the
data to a subdirectory named after the instrument, use

    >>> res=client.get(qr, path='/Users/ireland/Desktop/Data/{instrument}/{file}.fits')

Note that the download process is spawned in parallel to your existing
Python session.  This means that the remainder of your Python script
will continue as the download proceeds.  This may cause a problem if
the remainder of your script relies on the presence of the downloaded
data.  If you want to resume your script after all the data has been
downloaded then append '.wait()' to the 'get' command above, i.e.,

     >>> res=client.get(qr, path='/Users/ireland/Desktop/Data/{instrument}/{file}.fits').wait()

More information on the options available can be found through the
standard Python 'help' command.

Using the legacy query keywords it is very easy to translate a
Solarsoft/IDL VSO command into the equivalent SunPy VSO legacy query.
However, more powerful queries are possible with the new query style,
which is described below.


4. The new query style
----------------------

The new query style makes more complex queries possible.  Let's start
with translating the above legacy query into the syntax of the new
query:

    >>> qr=client.query(vso.attrs.Time('2001/1/1', '2001/1/2'), vso.attrs.Instrument('eit'))

Let's break down the arguments of client.query.  The first argument:

    vso.attrs.Time('2001/1/1', '2001/1/2')

sets the start and end times for the query (as with the legacy query,
any date/time format understood by SunPy's 'parse_time' function can
be used to specify dates and time).  The second argument:

    vso.attrs.Instrument('eit')

sets the instrument we are looking for.  So what is going on here?
The notion is that a VSO query has a set of attribute objects -
described in 'vso.attrs' - that are specifed to construct the query.
For the full list of vso attributes, use

    >>> help(vso.attrs)

Note that due to quirks at the VSO, we do not recommend that the
extent object 'vso.attrs.Extent' be in your query.  Instead, we
recommend that any extent filtering you need to do be done on the
queries made without setting a value to the vso.attrs.Extent object.
As we will see, the new-style query can take more than two arguments,
each argument separated from the other by a comma.  Each of those
arguments are chained together using a logical "AND".

The new-style query allows you to combine these VSO attribute objects
in complex ways that are not possible with the legacy query style.

So, let's look for the EIT and MDI data on the same day:

    >>> qr=client.query(vso.attrs.Time('2001/1/1', '2001/1/2'), vso.attrs.Instrument('eit') | vso.attrs.Instrument('mdi'))
    >>> qr.num_records()
    233
    >>> qr.show()

The two instrument types are joined together by the operator '|'.
This is the 'or' operator.  Think of the above query as setting a set
of conditions which get passed to the VSO.  Let's say you want all the
EIT data from two separate days:

    >>> qr=client.query(vso.attrs.Time('2001/1/1', '2001/1/2') | vso.attrs.Time('2007/8/9', '2007/8/10'), vso.attrs.Instrument('eit') )
    >>> qr.num_records()
    227

Each of the arguments in the new-style query can be thought of as
setting conditions that the returned records must satisfy.  You can
set the wavelength; for example, to return the 171 Angstrom EIT results

    >>> qr=client.query(vso.attrs.Time('2001/1/1', '2001/1/2'), vso.attrs.Instrument('eit'), vso.attrs.Wave(171,171) )
    >>> qr.num_records()
    4

The new-style query returns the same type of response as the legacy
query.  This means you can use the same command and syntax as shown
above to download your data.

Finally, please let us know if you encounter a bug while using the VSO
capabilities of SunPy.  Bugs are best reported via the issue tracking
system on GitHub - go to https://github.com/sunpy/sunpy/issues and
click on the "New Issue" button.
