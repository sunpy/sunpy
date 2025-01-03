.. _sunpy-topic-guide-deprecation-versioning:

**********************************************
Deprecation, Versioning and Release Practices
**********************************************

Versioning Policy
=================

Sunpy uses a date-based versioning system with three components, each serving a specific purpose:

- **X (LTS Version Number)**: Incremented with every Long Term Support (LTS) release.
- **Y (Release Counter)**: Starts at 0 for LTS releases and increases for intermediate releases.
- **Z (Bug Fix Counter)**: Incremented for bug fix releases.

Here's an example release sequence to illustrate:

- **1.0.0**: An LTS release (November 2024).
- **1.0.1**: A bug fix for the LTS release (December 2024).
- **1.1.0**: A short-support release (May 2025).
- **1.0.2**: A bug fix for the LTS release (June 2025).
- **1.1.1**: A bug fix for the short-support release (June 2025).
- **2.0.0**: An LTS release (November 2025).

For detailed information, refer to `SunPy Enhancement Proposal (SEP) 0009 <https://github.com/sunpy/sunpy-SEP/blob/main/SEP-0009.md#deprecations-and-documentation>`__.
Release Practices
=================

There will be two planned major releases of ``sunpy`` per year with at around 6 months between these releases.
Users can expect these releases occur in the months of June and December based on the current release schedule of the ``astropy``.
``SunPy`` provides two types of releases:

- **Long Term Support (LTS) Releases**: These are first major release of each year and supported for 12 months or until the next LTS release.
- **Short-Support (Non-LTS) Releases**: These are the second release of the year, supported for 6 months or until the next release.

Support periods may be extended beyond these requirements if necessary.

Python and Package Support Policies
===================================

SunPy aims to support a wide range of Python versions:

- SunPy supports Python versions released within the past three years.
- Support for newly released Python versions may take a few months to ensure compatibility.

Deprecation Process
===================

Sunpy's deprecation policy ensures users have ample time to adapt to changes and removals in the API. Below are the key practices:

Deprecation Policy and Breaking Changes
---------------------------------------

Sunpy enforces strict standards for updating, changing, or deprecating public APIs across the core package and stable affiliated packages. If functionality is planned for deprecation, the following steps are taken:

1. Use the `@deprecated` decorator to mark code for deprecation.
   Example:

   .. code-block:: python

       from sunpy.util.decorators import deprecated

       @deprecated(since="6.1", message="We will be moving this", alternative="sunpy.net.Scraper")
       class Scraper:
           pass

2. For changes to default values of keyword arguments, issue a deprecation warning.
   Example:

   .. code-block:: python

       from sunpy.util.exceptions import warn_deprecated

       if response_format is None:
           response_format = "legacy"
           warn_deprecated("The default response format from the VSO client will "
                           "be changing to 'table' in version 6.1. "
                           "To remove this warning set response_format='legacy' "
                           "to maintain the old behaviour or response_format='table' "
                           "to use the new behaviour.")

3. A deprecation warning must be introduced in one LTS release before the deprecated code is removed. For instance:
   - The warning is introduced in SunPy 6.1.
   - The code can only be removed in SunPy 7.1, after the 7.0 LTS release.

4. Add a **“deprecation”** entry to the changelog when introducing the warning and a **“removal”** entry when the code is removed.

5. Every effort is made to avoid breaking changes in public APIs. In cases where it is unavoidable, SunPy conducts code searches on GitHub or GitLab to assess the potential impact on openly available user code.

Documentation of Changes
------------------------

To help users adapt, all API changes are documented thoroughly:

- Deprecation warnings inform users of planned changes or removals, provided an alternative exists or the functionality is being removed entirely.
- All changes must be documented in the CHANGELOG and user documentation.
- When practical, provide a side-by-side comparison of old and new functionality.
- Include a **“breaking”** changelog entry with detailed guidance for updating user code.
