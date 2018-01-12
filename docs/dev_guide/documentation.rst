Documentation
=============

All code must be documented. Undocumented code will not be accepted into SunPy.
Documentation should follow the guidelines in `PEP 8
<https://www.python.org/dev/peps/pep-0008/>`_ and `PEP 257 (Docstring
conventions) <https://www.python.org/dev/peps/pep-0257/>`_. Documentation for
modules, classes, and functions should follow the `NumPy/SciPy documentation
style guide
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.
We provide an example of good documentation below or you can just browse some
of SunPy code itself for examples. All of the SunPy documentation
(like this page!) is built by Sphinx and must therefore adhere to Sphinx
guidelines.

Sphinx
------

`Sphinx <http://www.sphinx-doc.org/en/stable/>`_ is a tool for generating high-quality
documentation in various formats (HTML, pdf, etc) and is especially well-suited
for documenting Python projects. Sphinx works by parsing files written using a
`a Mediawiki-like syntax
<http://docutils.sourceforge.net/docs/user/rst/quickstart.html>`_ called
`reStructuredText <http://docutils.sourceforge.net/rst.html>`_. In addition
to parsing static files of reStructuredText, Sphinx can also be told to parse
code comments. In fact, in addition to what you are reading right now, the
`Python documentation <https://www.python.org/doc/>`_ was also created using
Sphinx.

Usage
#####

All of the SunPy documentation is contained in the ``doc/source`` folder and
code comments. The examples from the example gallery can be found in
``examples``. To build the documentation locally you must have Sphinx
(as well as Numpydoc, astropy-helpers, and sphinx-gallery) installed on
your computer. In the root directory run ::

    python setup.py build_docs

This will generate HTML documentation for SunPy in the ``docs/_build/html``
directory. The gallery examples are located under
``docs/_build/html/generated/gallery`` Sphinx builds documentation
iteratively only adding things that have changed. If you'd like to start
from scratch then just delete the build directory.

For more information on how to use Sphinx, consult the `Sphinx documentation
<http://www.sphinx-doc.org/en/stable/contents.html>`_.

The rest of this section will describe how to document the SunPy code in order
to guarantee well-formatted documentation.

doctest
#######

The example codes in the Guide section of the docs are configured with the Sphinx
`doctest extension <http://www.sphinx-doc.org/en/stable/ext/doctest.html>`_.
This will test the example code to make sure it runs correctly, it can be executed
using: ::

  sphinx-build -t doctest -b doctest ./ ../build

from inside the ``doc/source`` folder.

Use of quantities and units
---------------------------

Much code perform calculations using physical quantities.  SunPy uses astropy's
`quantities and units <http://docs.astropy.org/en/stable/units/index.html>`__
implementation to store, express and convert physical quantities. New classes
and functions should adhere to SunPy's `quantity and unit usage guidelines
<https://github.com/sunpy/sunpy-SEP/blob/master/SEP-0003.md>`__.  This document
sets out SunPy's reasons and requirements for the usage of quantities and
units.  Briefly, SunPy's `policy <https://github.com/sunpy/sunpy-SEP/blob/master/SEP-0003.md>`__
is that *all user-facing function/object arguments which accept physical
quantities as input **MUST** accept astropy quantities, and **ONLY** astropy
quantities*.

Developers should consult the
`Astropy Quantities and Units page <http://docs.astropy.org/en/stable/units/index.html>`__
for the latest updates on using quantities and units.  The `astropy tutorial on quantities and units
<http://www.astropy.org/astropy-tutorials/Quantities.html>`__ also provides useful examples on their
capabilities.

Astropy provides the decorator `~astropy.units.quantity_input` that
checks the units of the input arguments to a function against the
expected units of the argument.  We recommend using this decorator to
perform function argument unit checks.  The decorator ensures that the
units of the input to the function are convertible to that specified
by the decorator, for example ::

    >>> import astropy.units as u
    >>> @u.quantity_input(myangle=u.arcsec)
    ... def myfunction(myangle):
    ...     return myangle**2

This function only accepts arguments that are convertible to arcseconds.
Therefore, ::

    >>> myfunction(20 * u.degree)
    <Quantity 400. deg2>

returns the expected answer but ::

    >>> myfunction(20 * u.km)
    Traceback (most recent call last):
    ...
    astropy.units.core.UnitsError: Argument 'myangle' to function 'myfunction' must be in units convertible to 'arcsec'.

raises an error.

The following is an example of a use-facing function that returns the area of a
square, in units that are the square of the input length unit::

    >>> @u.quantity_input(side_length=u.m)
    ... def get_area_of_square(side_length):
    ...     """
    ...     Compute the area of a square.
    ...
    ...     Parameters
    ...     ----------
    ...     side_length : `~astropy.units.quantity.Quantity`
    ...         Side length of the square
    ...
    ...     Returns
    ...     -------
    ...     area : `~astropy.units.quantity.Quantity`
    ...         Area of the square.
    ...     """
    ...
    ...     return (side_length ** 2)

This more advanced example shows how a private function that does not accept
quantities can be wrapped by a function that does::

    >>> @u.quantity_input(side_length=u.m)
    ... def some_function(length):
    ...     """
    ...     Does something useful.
    ...
    ...     Parameters
    ...     ----------
    ...     length : `~astropy.units.quantity.Quantity`
    ...         A length.
    ...
    ...     Returns
    ...     -------
    ...     length : `~astropy.units.quantity.Quantity`
    ...         Another length
    ...     """
    ...
    ...     # the following function either
    ...     # a] does not accept Quantities
    ...     # b] is slow if using Quantities
    ...     result = _private_wrapper_function(length.convert('meters').value)
    ...
    ...     # now convert back to a quantity
    ...     result = Quantity(result_meters, units_of_the_private_wrapper_function)
    ...
    ...     return result

In this example, the non-user facing function *_private_wrapper_function* requires a numerical input in units of
meters, and returns a numerical output.  The developer knows that the result of *_private_wrapper_function* is in the
units *units_of_the_private_wrapper_function*, and sets the result of *some_function* to return the answer in those
units.

.. doctest-skip-all

Examples
--------

Modules
#######

Each module or package should begin with a docstring describing its overall
purpose and functioning. Below that meta-tags containing author, license, email
and credits information may also be listed.

Example: ::

    """This is an example module comment.

    An explanation of the purpose of the module would go here and will appear
    in the generated documentation
    """
    #
    # TODO
    #  Developer notes and todo items can be listed here and will not be
    #  included in the documentation.
    #
    __authors__ = ["Keith Hughitt", "Steven Christe", "Jack Ireland", "Alex Young"]
    __email__ = "keith.hughitt@nasa.gov"
    __license__ = "xxx"

For details about what sections can be included, see the section on `documenting
modules
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_ in the
NumPy/SciPy style guide.

Classes
#######

Class docstrings should include a clear and concise docstring explaining the
overall purpose of the class, required and optional input parameters, and the
return value. Additionally, notes, references and examples are encouraged.

Example (:class:`sunpy.map.Map`) ::

    """
    Map(data, header)

    A spatially-aware data array based on the SolarSoft Map object

    Parameters
    ----------
    data : numpy.ndarray, list
        A 2d list or ndarray containing the map data
    header : dict
        A dictionary of the original image header tags

    Attributes
    ----------
    header : dict
        A dictionary representation of the image header
    date : datetime
        Image observation time
    det : str
        Detector name
    inst : str
        Instrument name
    meas : str, int
        Measurement name. For AIA this is the wavelength of image
    obs : str
        Observatory name
    r_sun : float
        Radius of the sun
    name : str
        Nickname for the image type (e.g. "AIA 171")
    center : dict
        X and Y coordinate for the center of the sun in arcseconds
    scale: dict
        Image scale along the x and y axes in arcseconds/pixel

    Examples
    --------
    >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    >>> aia.T
    Map([[ 0.3125,  1.    , -1.1875, ..., -0.625 ,  0.5625,  0.5   ],
    [-0.0625,  0.1875,  0.375 , ...,  0.0625,  0.0625, -0.125 ],
    [-0.125 , -0.8125, -0.5   , ..., -0.3125,  0.5625,  0.4375],
    ...,
    [ 0.625 ,  0.625 , -0.125 , ...,  0.125 , -0.0625,  0.6875],
    [-0.625 , -0.625 , -0.625 , ...,  0.125 , -0.0625,  0.6875],
    [ 0.    ,  0.    , -1.1875, ...,  0.125 ,  0.    ,  0.6875]])
    >>> aia.header['cunit1']
    'arcsec'
    >>> aia.show()
    >>> import matplotlib.cm as cm
    >>> import matplotlib.colors as colors
    >>> aia.peek(cmap=cm.hot, norm=colors.Normalize(1, 2048))

    See Also
    --------
    numpy.ndarray Parent class for the Map object

    References
    ----------
    | http://docs.scipy.org/doc/numpy/reference/arrays.classes.html
    | http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    | https://www.scipy.org/Subclasses

    """

Functions
#########

Functions should include a clear and concise docstring explaining the overall
purpose of the function, required and optional input parameters, and the return
value. Additionally, notes, references and examples are encouraged.

Example (`numpy.matlib.ones
<https://github.com/numpy/numpy/blob/master/numpy/matlib.py>`_): ::

    def ones(shape, dtype=None, order='C'):
        """
        Matrix of ones.

        Return a matrix of given shape and type, filled with ones.

        Parameters
        ----------
        shape : {sequence of ints, int}
            Shape of the matrix
        dtype : data-type, optional
            The desired data-type for the matrix, default is np.float64.
        order : {'C', 'F'}, optional
            Whether to store matrix in C- or Fortran-contiguous order,
            default is 'C'.

        Returns
        -------
        out : matrix
            Matrix of ones of given shape, dtype, and order.

        See Also
        --------
        ones : Array of ones.
        matlib.zeros : Zero matrix.

        Notes
        -----
        If `shape` has length one i.e. ``(N,)``, or is a scalar ``N``,
        `out` becomes a single row matrix of shape ``(1,N)``.

        Examples
        --------
        >>> np.matlib.ones((2,3))
        matrix([[ 1.,  1.,  1.],
                [ 1.,  1.,  1.]])

        >>> np.matlib.ones(2)
        matrix([[ 1.,  1.]])

        """
        a = ndarray.__new__(matrix, shape, dtype, order=order)
        a.fill(1)
        return a

For details about what sections can be included, see the section on `documenting
functions
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_ in the
NumPy/SciPy style guide.

Trouble-shooting
----------------

Sphinx can be very particular about formatting, and the warnings and errors
outputted aren't always obvious.

Below are some commonly-encountered warning/error messages along with a
human-readable translation:

**WARNING: Duplicate explicit target name: "xxx".**

If you reference the same URL, etc more than once in the same document sphinx
will complain. To avoid, use double-underscores instead of single ones after
the URL.

**ERROR: Malformed table. Column span alignment problem at line offset n**

Make sure there is a space before and after each colon in your class and
function docs (e.g. attribute : type, instead of attribute: type). Also, for
some sections (e.g. Attributes) numpydoc seems to complain when a description
spans more than one line, particularly if it is the first attribute listed.

**WARNING: Block quote ends without a blank line; unexpected unindent.**

Lists should be indented one level from their parents.

**ERROR: Unkown target name: "xxx"**

In addition to legitimate errors of this type, this error will also occur when
variables have a trailing underscore, e.g., ``xxx_``.

**WARNING: Explicit markup ends without a blank line; unexpected unindent.**

This usually occurs when the text following a directive is wrapped to the next
line without properly indenting a multi-line text block.

**WARNING: toctree references unknown document '...'** /
**WARNING: toctree contains reference to nonexisting document**

This pair of errors is due to the way numpydoc scrapes class members.
