.. _docs_guidelines:

*************************
SunPy Documentation Rules
*************************

Overview
========

All code must be documented and we follow these style conventions described here:

* `numpydoc <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_
* `astropy <https://docs.astropy.org/en/latest/development/docrules.html>`_

We recommend familiarizing yourself with these references.

Differences
-----------

The current differences we have are as follows:

We backtick each type in the documentation strings so that they are interlinked by our documentation builder:

.. code-block:: python

    """
    Parameters
    ----------
    x : `type`
       Description of parameter x.
    """

SunPy Specific Rules
--------------------

* For **all** RST files, we enforce a one sentence per line rule and ignore the line length.

* Core datatypes are referenced as follows:

.. code-block:: python

    """
    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A SunPy map.
    """

* The Factory classes (``Map``, ``Timeseries``, ``Fido``) are referenced as follows:

.. code-block:: python

    `Map <sunpy.map.map_factory.MapFactory>`
    `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>`
    `TimeSeries <sunpy.timeseries.timeseries_factory.TimeSeriesFactory>`

* Standards on docstring length and style are enforced using `docformatter <https://pypi.org/project/docformatter/>`__:

.. code-block:: bash

    $ docformatter -r -i  --pre-summary-newline --make-summary-multi-line

.. _Docs Guidelines for Data Sources:

Documenting Data Sources
----------------------------
Subclasses of `~sunpy.map.GenericMap` or `~sunpy.timeseries.TimeSeries` must provide a detailed docstring providing an overview of the data source that the object represents.
In order to maintain consistency and completeness, the following information must be provided by a data source docstring, if available, and preferably in the following order:

* the name of the mission and instrument and the institution that built it
* short description of the instrument (e.g. Cassegrain reflector, Wolter-1 grazing incidence x-ray, coronagraph) including the type of detector
* description of the platform (e.g. satellite in 28 deg inclined orbit, a telescope on the summit of Mauna Kea in Hawaii)
* description of the primary purpose or science goals of the instrument.
* list of all wavelength(s) or passbands in appropriate units
* description of the emission processes which dominate in those passbands
* appropriate measurement properties such as field of view, angular resolution, time resolution
* description of the operational concept (e.g. operates 24/7, observes from 7 am to 5 pm UT) including mention of unusual operations scenarios (e.g. calibration seasons, eclipse seasons)
* the start and end of the data set

In addition, a reference section must be provided with links to the following resources, if available,

* the mission web page
* the instrument web page
* relevant wikipedia page(s)
* relevant user guide(s)
* the mission paper and instrument paper
* information to interpret metadata keywords such as FITS header reference
* the data archive

An example docstring can be found in the :ref:`Writing a new Instrument Map Class guide<new_maps_ts_etc>`.

Sphinx
======

All of the SunPy documentation (like this page) is built by `Sphinx <https://www.sphinx-doc.org/en/stable/>`_, which is a tool especially well-suited for documenting Python projects.
Sphinx works by parsing files written using a `a Mediawiki-like syntax <http://docutils.sourceforge.net/docs/user/rst/quickstart.html>`_ called `reStructuredText <http://docutils.sourceforge.net/rst.html>`_.
In addition to parsing static files of reStructuredText, Sphinx can also be told to parse code comments.
In fact, in addition to what you are reading right now, the `Python documentation <https://www.python.org/doc/>`_ was also created using Sphinx.

Usage
-----

All of the SunPy documentation is contained in the "docs" folder and code docstings/comments.
The examples from the example gallery can be found in the "examples" folder.
To build the documentation locally you must have all the dependencies (``pip install -e .[docs]``) specified in ``setup.cfg`` installed on your computer.

In the root directory run::

    $ tox -e build_docs

This will generate HTML documentation for SunPy in the "docs/_build/html" directory.
You can open the "index.html" file to browse the final product.
The gallery examples are located under "docs/_build/html/generated/gallery".
Sphinx builds documentation iteratively, only adding things that have changed.
If you'd like to start from scratch then just delete the tox build directory::

    $ rm .tox/build_docs/

or::

    $ tox --recreate -e build_docs

For more information on how to use Sphinx, consult the `Sphinx documentation <http://www.sphinx-doc.org/en/stable/contents.html>`_.

Trouble-shooting
----------------

Sphinx can be very particular about formatting, and the warnings and errors aren't always obvious.

Below are some commonly-encountered warning/error messages along with a human-readable translation:

**WARNING: Duplicate explicit target name: "xxx".**

If you reference the same URL, etc more than once in the same document sphinx will complain.
To avoid, use double-underscores instead of single ones after the URL.

**ERROR: Malformed table. Column span alignment problem at line offset n**

Make sure there is a space before and after each colon in your class and
function docs (e.g. attribute : type, instead of attribute: type).
Also, for some sections (e.g. Attributes) numpydoc seems to complain when a description spans more than one line, particularly if it is the first attribute listed.

**WARNING: Block quote ends without a blank line; unexpected unindent.**

Lists should be indented one level from their parents.

**ERROR: Unknown target name: "xxx"**

In addition to legitimate errors of this type, this error will also occur when variables have a trailing underscore, e.g., ``xxx_``.

**WARNING: Explicit markup ends without a blank line; unexpected unindent.**

This usually occurs when the text following a directive is wrapped to the next line without properly indenting a multi-line text block.

**WARNING: toctree references unknown document '...'** / **WARNING: toctree contains reference to nonexisting document**

This pair of errors is due to the way numpydoc scrapes class members.
