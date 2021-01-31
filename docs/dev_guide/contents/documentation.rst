.. _docs_guidelines:

*************************
SunPy Documentation Rules
*************************

Overview
========

All code must be documented and we follow these style conventions described here:

* `numpydoc <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_

We recommend familiarizing yourself with this style.

Referring to other code
-----------------------

To link to other methods, classes, or modules in sunpy you have to use backticks, for example:

.. code-block:: rst

    `sunpy.io.read_file`

generates a link like this: `sunpy.io.read_file`.

We use the sphinx setting ``default_role = 'obj'`` so that you do not nor **SHOULD NOT** use the ``:class:`` qualifier, but ``:func:``, ``:meth:`` are different (more on this below).

Often, you don't want to show the full package and module name.
As long as the target is unambiguous you can simply leave them out:

.. code-block:: rst

    `.read_file`

and the link still works: `.read_file`.

If there are multiple code elements with the same name (e.g. ``peek()`` is a method in multiple classes), you'll have to extend the definition:

.. code-block:: rst

    `.GenericMap.peek` or `.CompositeMap.peek`

These will show up as `.GenericMap.peek` or `.CompositeMap.peek`.
To still show only the last segment you can add a tilde as prefix:

.. code-block:: rst

    `~.GenericMap.peek` or `~.CompositeMap.peek`

will render as `~.GenericMap.peek` or `~.CompositeMap.peek`.

Other packages can also be linked via
`intersphinx <http://www.sphinx-doc.org/en/master/ext/intersphinx.html>`_:

.. code-block:: rst

    `numpy.mean`

will return this link: `numpy.mean`.
This works for Python, Numpy and Astropy (full list is in :file:`docs/conf.py`).

With Sphinx, if you use ``:func:`` or ``:meth:``, it will add closing brackets to the link.
If you get the wrong pre-qualifier, it will break the link, so we suggest that you double check if what you are linking is a method or a function.

.. code-block:: rst

    :class:`numpy.mean()`
    :meth:`numpy.mean()`
    :func:`numpy.mean()`

will return two broken links ("class" and "meth") but "func" will work:

"class": :class:`numpy.mean()`

"meth" :meth:`numpy.mean()`

"func" :func:`numpy.mean()`

SunPy-Specific Rules
--------------------

* For **all** RST files, we enforce a one sentence per line rule and ignore the line length.
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

An example docstring can be found in the :ref:`Writing a new Instrument Map Class guide <new_maps_ts_etc>`.

Sphinx
======

All of the SunPy documentation (like this page) is built by `Sphinx <https://www.sphinx-doc.org/en/stable/>`_, which is a tool especially well-suited for documenting Python projects.
Sphinx works by parsing files written using a `a Mediawiki-like syntax <http://docutils.sourceforge.net/docs/user/rst/quickstart.html>`_ called `reStructuredText <http://docutils.sourceforge.net/rst.html>`_.
In addition to parsing static files of reStructuredText, Sphinx can also be told to parse code comments.
In fact, in addition to what you are reading right now, the `Python documentation <https://www.python.org/doc/>`_ was also created using Sphinx.

Usage
-----

All of the SunPy documentation is contained in the "docs" folder and code documentation strings.
The examples from the example gallery can be found in the "examples" folder.

In the root directory run::

    $ tox -e build_docs

This will generate HTML documentation for SunPy in the "docs/_build/html" directory.
You can open the "index.html" file to browse the final product.
The gallery examples are located under "docs/_build/html/generated/gallery".
Sphinx builds documentation iteratively, only adding things that have changed.

If you want to build the documentation without building the gallery, i.e. to reduce build times while working on other sections of the documentation you can run::

    $ tox -e build_docs -- -D plot_gallery=False

If you'd like to start from scratch (i.e., remove the tox cache) then run::

    $ tox -e build_docs -- -aE

To build the documentation in your current python environment you must have all the dependencies specified in ``setup.cfg`` installed (``pip install -e .[docs]``).
Then change to the :file:`docs/` directory and run::

    $ make html

For more information on how to use Sphinx, consult the `Sphinx documentation <http://www.sphinx-doc.org/en/stable/contents.html>`_.

Special Sphinx directives
-------------------------

``minigallery`` directive
^^^^^^^^^^^^^^^^^^^^^^^^^

Sphinx will automatically record which functions, classes, etc. are used in each gallery example.
In the documentation, you can insert a mini-gallery of the subset of the gallery examples that uses a particular function, class, etc.
For example, the following RST block::

    .. minigallery:: sunpy.coordinates.RotatedSunFrame

produces this mini-gallery:

.. minigallery:: sunpy.coordinates.RotatedSunFrame

If you want to specify more than one object, separate them by spaces.
This is particularly useful if you need to cover multiple namespaces in which an object may be accessed, e.g.::

    .. minigallery:: sunpy.coordinates.RotatedSunFrame sunpy.coordinates.metaframes.RotatedSunFrame

``generate`` directive
^^^^^^^^^^^^^^^^^^^^^^

In rare circumstances, one may want to insert "raw" HTML directly into the pages written by Sphinx.
For HTML that is statically available (i.e., already written in some form), one can use the `"raw" directive <https://docutils.sourceforge.io/docs/ref/rst/directives.html#raw-data-pass-through>`__.
For HTML that is generated by Python code, SunPy provides the custom directive ``generate``.
Here's an example RST block::

    .. generate:: html
        :html_border:

        import os
        from sunpy.data.sample import file_dict
        print("<table>")
        for key, value in file_dict.items():
            print(f"<tr><th>{key}</th><td>{os.path.basename(value)}</td></tr>")
        print("</table>")

to insert the following HTML table:

.. generate:: html
    :html_border:

    import os
    from sunpy.data.sample import file_dict
    print("<table>")
    for key, value in file_dict.items():
        print(f"<tr><th>{key}</th><td>{os.path.basename(value)}</td></tr>")
    print("</table>")

Troubleshooting
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
