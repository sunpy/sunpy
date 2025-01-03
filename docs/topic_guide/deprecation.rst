.. _sunpy-topic-guide-logging-system:

********************************************
Deprecation Versioning and Release Practices
********************************************

Versioning Policy
=================

SunPy follows a date-based versioning system. For detailed information, refer to `SunPy Enhancement Proposal (SEP) 0009 <https://github.com/sunpy/sunpy-SEP/blob/main/SEP-0009.md#deprecations-and-documentation>`_.

Release Practices
=================

We aim to have a few minor releases per year, though the exact number may vary based on available funding and resources.

Python and package support policies
=====================================
Support versions of Python released in last three years
Support for newest Python may take a few months

Deprecation Process
===================

Our deprecation policy is designed to help users prepare for changes and removals in SunPy. This section explains why you might see deprecation warnings, what they mean, and how you can address them.

Deprecation Policy and Breaking Changes
---------------------------------------

All public API within the Sunpy project (the sunpy core package and stable affiliated packages) will enforce strict standards when it comes to either changing, updating or breaking the API.
When we plan to deprecate functionality, the following steps are taken to ensure users have ample time to adapt:

If you want to deprecate anything within in sunpy, you should do the following:

.. code-block:: python

    from sunpy.util.decorators import deprecated

    @deprecated(since="3.1", message="We will be moving this", alternative="sunpy.net.Scraper")
    class Scraper:


If changing the default value of a keyword argument, a warning is issued.
Example:

.. code-block:: python

    from sunpy.util.exceptions import warn_deprecated

    if response_format is None:
        response_format = "legacy"
        warn_deprecated("The default response format from the VSO client will "
                        "be changing to 'table' in version 3.1. "
                        "To remove this warning set response_format='legacy' "
                        "to maintain the old behaviour or response_format='table'"
                        " to use the new behaviour.")

The deprecation warning has to be in one LTS release before the deprecated code can be removed. So in the above example, the warning will be in sunpy 3.1 but it can not be removed until sunpy 4.1 after the 4.0 LTS release.

A “deprecation” entry is added to the changelog when the warning is introduced. A “removal” entry is added when the deprecated code is finally removed.

Documentation of Changes
------------------------
Every attempt is made to avoid breaking any public API in sunpy but in the case it does happen.

There should be a “breaking” changelog entry to accompany the change with as much detail on how to update a user’s code due to this breaking change.

If user code is made openly available on GitHub or GitLab, we can do code searches to see if any deprecations or removals would have a major impact.


