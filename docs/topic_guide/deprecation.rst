.. _sunpy-topic-guide-deprecation-versioning:

**********************************************
Deprecation, Versioning and Release Practices
**********************************************

Versioning Policy
=================

Sunpy uses semantic versioning system with three components, each serving a specific purpose:

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

For detailed information, refer to `Sunpy Enhancement Proposal (SEP) 0009 <https://github.com/sunpy/sunpy-SEP/blob/main/SEP-0009.md#deprecations-and-documentation>`__.

Release Practices
=================

There will be two planned major releases of ``sunpy`` per year, with approximately 6 months between these releases. Users can expect these releases to occur in the months of June and December, based on the current release schedule of ``astropy``.

``sunpy`` provides two types of releases:

- **Long Term Support (LTS) Releases**: These are the first major release of each year and are supported for 12 months or until the next LTS release.
- **Short-Support (Non-LTS) Releases**: These are the second release of the year, supported for 6 months or until the next release.

Support periods may be extended beyond these requirements if necessary.

Python and Package Support Policies
===================================

Sunpy aims to support a wide range of Python versions:

- Sunpy supports Python versions released within the past three years.
- Support for newly released Python versions may take a few months to ensure compatibility.

Sunpy has a short list of core dependencies (Python, numpy, astropy, parfive) and a long list of optional dependencies. The minimum version of these packages that we enforce follows this policy:

- **Python**: Released in the prior 36 months from the anticipated release date.
- **astropy**: Released in the prior 12 months from the anticipated release date.
- **Everything else**: Released in the prior 24 months from the anticipated release date.

Affiliated packages maintained by the Sunpy team will follow this policy and will support at least the Sunpy LTS version at the time of their release.

For dependencies only needed to run tests, Sunpy will support versions released in the prior 12 months to the current date, unless there is a critical issue that requires a newer version.

Deprecation Process
===================

Sunpy's deprecation policy ensures users have ample time to adapt to changes and removals in the API. Below are the key practices:

Deprecation Policy and Breaking Changes
---------------------------------------

Sunpy enforces strict standards for updating, changing, or deprecating public APIs across the core package and stable affiliated packages. If functionality is planned for deprecation, the following steps are taken:

- A deprecation warning must be introduced in one LTS release before the deprecated code is removed. For instance:
  - The warning is introduced in Sunpy 6.1.
  - The code can only be removed in Sunpy 7.1, after the 7.0 LTS release.

- Deprecation warnings inform users of planned changes or removals, provided an alternative exists or the functionality is being removed entirely.
- Include a **“breaking”** changelog entry with detailed guidance for updating user code.

Every effort is made to avoid breaking changes in public APIs. In cases where it is unavoidable, Sunpy conducts code searches on GitHub or GitLab to assess the potential impact on openly available user code.

When practical, provide a side-by-side comparison of old and new functionality.
