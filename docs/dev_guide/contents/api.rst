.. _public_api:

******************
sunpy's Public API
******************

Convention in the Python ecosystem is to add an underscore to the start of a function to denote if a function or method is "private" e.g., `~sunpy.coordinates.sun._angular_radius`.
If it is considered to be private, there are no guarantee that the API will change with a warning and so external use of these functions or methods are strongly discouraged.

sunpy follows this convention but with one extra caveat.
Within each python file, we have a ``__all__`` that defines what is imported into the namespace if you do e.g., ``from sunpy.coordinates.sun import *``.
This is the "public" API of that module or part of sunpy as well.
These functions are listed within our API documentation: :ref:`reference`.
If you do ``import sunpy.coordinates.sun``, you can still access the "private" functions.

This means that all of this code will follow the deprecation policy detailed below with the exception of `sunpy.util` which is considered for internal sunpy use only.

Deprecation Policy and Breaking Changes
=======================================

All public API within the SunPy project (sunpy and its affiliated packages) will enforce strict standards when it comes to either changing, updating or breaking the API.

.. _deprecation:

Deprecations
------------

If you want to deprecate anything within in sunpy, you should do the following:

.. code-block:: python

    from sunpy.util.decorators import deprecated

    @deprecated(since="3.1", message="We will be moving this", alternative="sunpy.net.Scraper")
    class Scraper:

The deprecation warning has to be in one LTS release before the deprecated code can be removed.
So in the above example, the warning will be in sunpy 3.1 but it can not be removed until sunpy 4.1 after the 4.0 LTS release.

There should be a "deprecation" changelog entry to accompany the deprecation warning.
When the code is actually removed, a "removal" changelog will be added.

The same applies if you want to change the default of an argument or keyword for a function or method.

.. code-block:: python

    if response_format is None:
        response_format = "legacy"
        warnings.warn("The default response format from the VSO client will "
                    "be changing to 'table' in version 3.1. "
                    "To remove this warning set response_format='legacy' "
                    "to maintain the old behaviour or response_format='table'"
                    " to use the new behaviour.",
                    SunpyDeprecationWarning,
                    stacklevel=2)

.. note::

    This is a summary of `SEP-0009`_ which as the formal SunPy deprecation policy.

.. _SEP-0009: https://github.com/sunpy/sunpy-SEP/blob/master/SEP-0009.md#deprecations-and-documentation

.. _breaking:

Breaking Changes
----------------

Every attempt is made to avoid breaking any public API in sunpy but in the case it does happen.

There should be a "breaking" changelog entry to accompany the change with as much detail on how to update a user's code due to this breaking change.
