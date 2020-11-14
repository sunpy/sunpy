.. doctest-skip-all

.. _units_in_code:

***************************
Use of quantities and units
***************************

Much code perform calculations using physical quantities.
SunPy uses astropy's `quantities and units <https://docs.astropy.org/en/stable/units/index.html>`_ implementation to store, express and convert physical quantities.
New classes and functions should adhere to SunPy's `quantity and unit usage guidelines
<https://github.com/sunpy/sunpy-SEP/blob/master/SEP-0003.md>`_.

This document sets out SunPy's reasons and requirements for the usage of quantities and units.
Briefly, SunPy's `policy <https://github.com/sunpy/sunpy-SEP/blob/master/SEP-0003.md>`_ is that *all user-facing function/object arguments which accept physical quantities as input **MUST** accept astropy quantities, and **ONLY** astropy quantities*.

Developers should consult the `Astropy Quantities and Units page <https://docs.astropy.org/en/stable/units/index.html>`_ for the latest updates on using quantities and units.  The `astropy tutorial on quantities and units <https://www.astropy.org/astropy-tutorials/Quantities.html>`_ also provides useful examples on their
capabilities.

Astropy provides the decorator `~astropy.units.quantity_input` that checks the units of the input arguments to a function against the expected units of the argument.
We recommend using this decorator to perform function argument unit checks.
The decorator ensures that the units of the input to the function are convertible to that specified by the decorator, for example ::

    >>> import astropy.units as u
    >>> @u.quantity_input
    ... def myfunction(myangle: u.arcsec):
    ...     return myangle**2

This function only accepts arguments that are convertible to arcseconds.
Therefore::

    >>> myfunction(20 * u.degree)
    <Quantity 400. deg2>

returns the expected answer but::

    >>> myfunction(20 * u.km)
    Traceback (most recent call last):
    ...
    astropy.units.core.UnitsError: Argument 'myangle' to function 'myfunction' must be in units convertible to 'arcsec'.

raises an error.

The following is an example of a use-facing function that returns the area of a square, in units that are the square of the input length unit::

    >>> @u.quantity_input
    ... def get_area_of_square(side_length: u.m):
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

This more advanced example shows how a private function that does not accept quantities can be wrapped by a function that does::

    >>> @u.quantity_input
    ... def some_function(length: u.m):
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

In this example, the non-user facing function ``_private_wrapper_function`` requires a numerical input in units of meters, and returns a numerical output.
The developer knows that the result of ``_private_wrapper_function`` is in the units ``units_of_the_private_wrapper_function``, and sets the result of ``some_function`` to return the answer in those units.
