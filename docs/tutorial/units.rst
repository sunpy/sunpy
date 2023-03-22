.. _units-sunpy:

*****
Units
*****

This section of the guide talks about representing physical units in sunpy.
sunpy makes use of the `astropy.units` for this task.

All functions in sunpy that accept or return numbers associated with physical quantities accept and return `~astropy.units.Quantity` objects.
These objects represent a number (or an array of numbers) and a unit.
This means sunpy is always explicit about the units associated with a value.
Quantities and units are powerful tools for keeping track of variables with a physical meaning and make it straightforward to convert the same physical quantity into different units.

To use units we must first import them from Astropy.
To save on typing it's standard practice to import the units module as ``u``::

   >>> import astropy.units as u

Adding Units to Data
====================

We can create a `~astropy.units.Quantity` by multiplying a number by a unit::

   >>> length = 10 * u.meter
   >>> length
   <Quantity 10. m>

A `~astropy.units.Quantity` has both a ``.unit`` and a ``.value`` attribute::

  >>> length.value
  10.0

  >>> length.unit
  Unit("m")

Arithmetic Operations With Units
================================

Probably most usefully, `~astropy.units.Quantity` objects will propagate units through arithmetic operations when appropriate::

  >>> distance_start = 10 * u.mm
  >>> distance_end = 23 * u.km
  >>> length = distance_end - distance_start
  >>> length
  <Quantity 22.99999 km>

  >>> time = 15 * u.minute
  >>> speed = length / time
  >>> speed
  <Quantity 1.53333267 km / min>

However, operations which do not make physical sense for the units specified will cause an error::

  >>> length + time
  Traceback (most recent call last):
  ...
  astropy.units.core.UnitConversionError: Can only apply 'add' function to quantities with compatible dimensions

Converting Units
================

These `~astropy.units.Quantity` objects can also be converted to other units or unit systems::

  >>> length.to(u.km)
  <Quantity 0.01 km>

  >>> length.cgs
  <Quantity 1000. cm>

Unit Equivalencies
==================

Dropping Units
==============

As shown above, you can retrieve the just the numerical value of a `~astropy.units.Quantity`,

  >>> length.value


Quantities as function arguments
================================

An extremely useful addition to the base functionality of Quanitities is the ``@u.quantity_input`` decorator.
This allows specification for required units for function arguments to ensure that the calculation within that function always makes physical sense.
For instance, if we defined a function to calculate speed as above, we might want the distance and time as inputs::

  >>> def speed(length, time):
  ...     return length / time

However, this assumes that the length and time passed in always have the appropriate units.
To enforce the correct units we can use `~astropy.units.quantity_input`, with `function annotations <https://python-3-for-scientists.readthedocs.io/en/latest/python3_features.html#function-annotations>`__ to specify the units::

  >>> @u.quantity_input
  ... def speed(length: u.m, time: u.s):
  ...     return length / time

Now, when this function is called, if the units of length and time are not convertible to the units specified, an error will be raised stating that the units are incorrect or missing::

  >>> speed(1*u.m, 10*u.m)
  Traceback (most recent call last):
  ...
  astropy.units.core.UnitsError: Argument 'time' to function 'speed' must be in units convertible to 's'.

  >>> speed(1*u.m, 10)
  ...
  Traceback (most recent call last):
  ...
  TypeError: Argument 'time' to function 'speed' has no 'unit' attribute. ... pass in an astropy Quantity instead.

Note that the units of the inputs do not have to be exactly the same as those in the function definition, as long as they can be converted to those units.
So for instance, passing in a time in minutes still works even though we specified ``time: u.s``::

  >>> speed(1*u.m, 1*u.minute)
  <Quantity 1. m / min>

This may still not be quite as we want it, since we wanted the input time in seconds but the output is in m/min.
We can correct this by defining the function with an additional annotation::

  >>> @u.quantity_input
  ... def speed(length: u.m, time: u.s) -> u.m/u.s:
  ...     return length / time

This will force the output of the function to be converted to m/s before returning, so that you will always have the same units on the output from this function::

  >>> speed(1*u.m, 1*u.minute)
  <Quantity 0.01666667 m / s>
