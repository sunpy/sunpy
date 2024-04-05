.. _sunpy-tutorial-units:

*****
Units
*****

In this section of the tutorial you will learn about representing physical units in sunpy.
All functions in sunpy that accept or return numbers associated with physical quantities do so using `astropy.units.Quantity` objects.
These objects represent a number (or an array of numbers) and a unit.
This means sunpy is always explicit about the units associated with a value.
Quantities and units are powerful tools for keeping track of variables with a physical meaning and make it straightforward to convert the same physical quantity into different units.

By the end of this section of the tutorial, you will learn how to create a `~astropy.units.Quantity`, perform basic arithmetic with a `~astropy.units.Quantity`, convert a `~astropy.units.Quantity` to different units, and write a function that ensures the inputs have the correct units.

Adding Units to Data
====================

To use units we must first import them from Astropy.
It is standard practice to import the units module as ``u``:

.. code-block:: python

   >>> import astropy.units as u

We can create a `~astropy.units.Quantity` by multiplying a number by a unit:

.. code-block:: python

   >>> length = 10 * u.meter
   >>> length
   <Quantity 10. m>

A `~astropy.units.Quantity` can be decomposed into its unit and numerical value using the ``.unit`` and ``.value`` attributes:

.. code-block:: python

   >>> length.value  # doctest: +SKIP
   10.0

   >>> length.unit
   Unit("m")

Arithmetic With Units
=====================

`~astropy.units.Quantity` objects propagate units through arithmetic operations:

.. code-block:: python

   >>> distance_start = 10 * u.mm
   >>> distance_end = 23 * u.km
   >>> displacement = distance_end - distance_start
   >>> displacement
   <Quantity 22.99999 km>

   >>> time = 15 * u.minute
   >>> speed = displacement / time
   >>> speed
   <Quantity 1.53333267 km / min>

However, operations with incompatible units raise an error:

.. code-block:: python

   >>> displacement + time
   Traceback (most recent call last):
   ...
   astropy.units.core.UnitConversionError: Can only apply 'add' function to quantities with compatible dimensions

Converting Units
================

`~astropy.units.Quantity` objects can also be converted to other units or unit systems:

.. code-block:: python

   >>> length.to(u.km)
   <Quantity 0.01 km>

   >>> length.cgs
   <Quantity 1000. cm>

Unit Equivalencies
==================

It is commonplace to convert between units which are only compatible under certain assumptions.
For example, in spectroscopy, spectral energy and wavelength are equivalent given the relation :math:`E=hc/\lambda`.
If we try to convert a wavelength to energy using what we learned in the previous section, we get an exception because length and energy are, in general, not compatible units:

.. code-block:: python

   >>> length.to(u.keV)
   Traceback (most recent call last):
   ...
   astropy.units.core.UnitConversionError: 'm' (length) and 'keV' (energy/torque/work) are not convertible

However, we can perform this conversion using the `~astropy.units.equivalencies.spectral` equivalency:

.. code-block:: python

   >>> length.to(u.keV, equivalencies=u.spectral())
   <Quantity 1.23984198e-10 keV>

An equivalency common in solar physics is conversion of angular distances in the plane of the sky to physical distances on the Sun.
To perform this conversion, sunpy provides `~sunpy.coordinates.utils.solar_angle_equivalency`, which requires specifying the location at which that angular distance was measured:

.. code-block:: python

   >>> from sunpy.coordinates import get_earth
   >>> from sunpy.coordinates.utils import solar_angle_equivalency

   >>> length.to(u.arcsec, equivalencies=solar_angle_equivalency(get_earth("2013-10-28")))
   INFO: Apparent body location accounts for 495.82 seconds of light travel time [sunpy.coordinates.ephemeris]
   <Quantity 1.38763748e-05 arcsec>

Note that in the above example we made use of `sunpy.coordinates.get_earth`.
We will talk more about coordinates in the :ref:`sunpy-tutorial-coordinates` section of this tutorial.
For now, it is just important to know that this function returns the location of the Earth on 2013 October 28.

Dropping Units
==============

Not every package in the scientific Python ecosystem understands units.
As such, it is sometimes necessary to drop the units before passing `~astropy.units.Quantity` to such functions.
As shown above, you can retrieve the just the numerical value of a `~astropy.units.Quantity`:

.. code-block:: python

   >>> length.to_value()  # doctest: +SKIP
   10.0
   >>> length.to_value(u.km)  # doctest: +SKIP
   0.01

Quantities as function arguments
================================

When calling a function that relies on inputs corresponding to physical quantities, there is often an implicit assumption that these input arguments are expressed in the expected units of that function.
For instance, if we define a function to calculate speed as above, the inputs should correspond to a distance and a time:

.. code-block:: python

   >>> def speed(length, time):
   ...     return length / time

However, this assumes that the two arguments passed in have units consistent with distance and time without checking.
The `~astropy.units.quantity_input` decorator, combined with `function annotations <https://python-3-for-scientists.readthedocs.io/en/latest/python3_features.html#function-annotations>`__, enforces compatible units on the function inputs:

.. code-block:: python

   >>> @u.quantity_input
   ... def speed(length: u.m, time: u.s):
   ...     return length / time

Now when this function is called, if the inputs are not convertible to the units specified, an error will be raised stating that the units are incorrect or missing:

.. code-block:: python

   >>> speed(1*u.m, 10*u.m)
   Traceback (most recent call last):
   ...
   astropy.units.core.UnitsError: Argument 'time' to function 'speed' must be in units convertible to 's'.

   >>> speed(1*u.m, 10)
   ...
   Traceback (most recent call last):
   ...
   TypeError: Argument 'time' to function 'speed' has no 'unit' attribute. ... pass in an astropy Quantity instead.

The units of the inputs need only be compatible with those in the function definition.
For example, passing in a time in minutes still works even though we specified ``time: u.s``:

.. code-block:: python

   >>> speed(1*u.m, 1*u.minute)
   <Quantity 1. m / min>

Note that the units of the output are dependent on the units of the inputs.
To ensure consistent units on the output of our function, we add an additional function annotation to force the output to always be converted to m/s before returning an answer:

.. code-block:: python

   >>> @u.quantity_input
   ... def speed(length: u.m, time: u.s) -> u.m/u.s:
   ...     return length / time
   >>> speed(1*u.m, 1*u.minute)
   <Quantity 0.01666667 m / s>
