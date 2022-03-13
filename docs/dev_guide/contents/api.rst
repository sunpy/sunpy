.. _public_api:

******************
sunpy's Public API
******************

Convention in the Python ecosystem is to add an underscore to the start of a function to denote if a function or method is "private" e.g., `~sunpy.coordinates.sun._angular_radius`.
If it is considered to be private, there is no guarantee that the API or behavior will change with a warning, so external use of these functions or methods are strongly discouraged.

sunpy follows this convention but with one extra caveat.
Within each python file, we have a ``__all__`` that defines what is imported into the namespace if you do e.g., ``from sunpy.coordinates.sun import *``.
This is the "public" API of that module.
These functions are the ones listed within our API documentation: :ref:`reference`.
If you do ``import sunpy.coordinates.sun``, you can still access the "private" functions.

This means that all of the public API will follow the deprecation policy detailed below with the exception of `sunpy.util` which is considered to be for internal sunpy use only.

Deprecation Policy and Breaking Changes
=======================================

All public API within the SunPy project (the sunpy core package and stable affiliated packages) will enforce strict standards when it comes to either changing, updating or breaking the API.

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
When the code is actually removed, a "removal" changelog should be added.

The same applies if you want to change the default value of a keyword argument for a function or method, e.g.:

.. code-block:: python

    from sunpy.util.exceptions import warn_deprecated

    if response_format is None:
        response_format = "legacy"
        warn_deprecated("The default response format from the VSO client will "
                        "be changing to 'table' in version 3.1. "
                        "To remove this warning set response_format='legacy' "
                        "to maintain the old behaviour or response_format='table'"
                        " to use the new behaviour.")

.. note::

    This is a summary of `SEP-0009`_ which is the formal SunPy project deprecation policy.

.. _SEP-0009: https://github.com/sunpy/sunpy-SEP/blob/master/SEP-0009.md#deprecations-and-documentation

.. _breaking:

Breaking Changes
----------------

Every attempt is made to avoid breaking any public API in sunpy but in the case it does happen.

There should be a "breaking" changelog entry to accompany the change with as much detail on how to update a user's code due to this breaking change.
