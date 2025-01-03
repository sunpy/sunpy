
==================

*****************
Deprecation Policy
*****************

Overview
--------

Our deprecation policy is designed to help users prepare for changes and removals in SunPy. This section explains why you might see deprecation warnings, what they mean, and how you can address them.

Versioning Policy
-----------------

**Date-based Versioning**

SunPy follows a date-based versioning system. For detailed information, refer to `SunPy Enhancement Proposal (SEP) 0009 <https://docs.sunpy.org/en/latest/dev_guide/sep/sep-0009.html>`_.

**Release Practices**

We aim to have a few minor releases per year, though the exact number may vary based on available funding and resources.


Deprecation Process
-------------------

When we plan to deprecate functionality, the following steps are taken to ensure users have ample time to adapt:

**Deprecation Warnings**

The deprecation warning has to be in one LTS release before the deprecated code can be removed. So in the above example, the warning will be in sunpy 3.1 but it can not be removed until sunpy 4.1 after the 4.0 LTS release.
Example:

.. code-block:: python

    from sunpy.util.decorators import deprecated

    @deprecated(since="3.1", message="We will be moving this", alternative="sunpy.net.Scraper")
    class Scraper:

**Changelog Entries**

A “deprecation” entry is added to the changelog when the warning is introduced. A “removal” entry is added when the deprecated code is finally removed.

**Deprecating Default Values**

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

Documentation of Changes
------------------------

**API Changes**

All API changes must be documented in the CHANGELOG and in the user documentation.

**Side-by-side Comparisons**

Where practical, a side-by-side comparison of old and new functionality will be provided.

Emission of Warnings
--------------------

**Deprecation Warnings**

Code will emit deprecation warnings to inform users of planned changes or removals.

**Availability of Alternatives**

Warnings will only be emitted when an alternative option exists for the functionality, or the functionality is to be completely removed.
