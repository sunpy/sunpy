.. _public_api:

************************
What is Public in sunpy?
************************

Convention in the Python ecosystem is to add an underscore to the start of a function to denote if a function or method is "private" e.g., `~sunpy.coordinates.sun._angular_radius`.
If it is considered to be private, there are no guarantee that the API will change with a warning and so external use of these functions or methods are strongly discouraged.

sunpy follows this convention but with one extra caveat.
Within each python file, we have a ``__all__`` that defines what is imported into the namespace if you do e.g., ``from sunpy.coordinates.sun import *``.
This is the "public" API of that module or part of sunpy as well.
These functions are listed within our API documentation: :ref:`reference`.
If you do ``import sunpy.coordinates.sun``, you can still access the "private" functions.

The only special case to this is `sunpy.util`, which while has a "public" API is not for public use and is used internally within sunpy and its affiliated packages.
