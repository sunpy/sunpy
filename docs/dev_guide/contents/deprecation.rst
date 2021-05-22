.. _deprecation:

******************
Deprecation Policy
******************

.. note::

    This page is a summary of `SEP-0009`_ which as the formal SunPy deprecation policy.


If you want to deprecate anything within in sunpy, you should do the following:

.. code-block:: python

    from sunpy.util.decorators import deprecated

    @deprecated(since="3.1", message="We will be moving this", alternative="sunpy.net.Scraper")
    class Scraper:

The deprecation warning has to be in one LTS release before the deprecated code can be removed.
So in the above example, the warning will be in sunpy 3.1 but it can not be removed until sunpy 4.1 after the 4.0 LTS release.

There should be a "breaking" changelog entry to accompany the deprecation warning.
When the code is actually removed, a "removal" changelog will be added.

.. _SEP-0009: https://github.com/sunpy/sunpy-SEP/blob/master/SEP-0009.md#deprecations-and-documentation
