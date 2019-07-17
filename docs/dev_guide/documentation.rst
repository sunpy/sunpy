.. _docs_guidelines:

*************************
SunPy Documentation Rules
*************************

Overview
========

All code must be documented and we follow these style conventions described here:

* `numpydoc <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_
* `astropy <https://docs.astropy.org/en/latest/development/docrules.html>`_

We recommend having a read of these and unless stated overwise, we stick to these standards.

Differences
-----------

The current differences we have is as follows:

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

* For all rst files, we enforce a one sentenace per line rule and ignore the line length.

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

    $ python setup.py build_docs

This will generate HTML documentation for SunPy in the "docs/_build/html" directory.
You can open the "index.html" file to browse the final product.
The gallery examples are located under "docs/_build/html/generated/gallery".
Sphinx builds documentation iteratively, only adding things that have changed.
If you'd like to start from scratch then just delete the build directory or run::

    $ python setup.py build_docs -l

to clean previous builds before building new ones.

For more information on how to use Sphinx, consult the `Sphinx documentation <http://www.sphinx-doc.org/en/stable/contents.html>`_.

Trouble-shooting
----------------

Sphinx can be very particular about formatting, and the warnings and errors outputted aren't always obvious.

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
