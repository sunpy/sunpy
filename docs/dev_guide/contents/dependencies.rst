.. _dependency_versions:

***********************
Dependencies and Extras
***********************

Dependency Support Policy
=========================

.. note::

    This policy is based on `NEP-0029`_.

sunpy has a short list of core dependencies (Python, numpy, astropy, parfive) and a long list of optional dependencies.
The minimum version of these packages that we enforce follows this policy.

* Python: Released in the prior 42 months from the anticipated release date.
* astropy: Released in the prior 12 months from the anticipated release date.
* Everything else: Released in the prior 24 months from the anticipated release date.

Sponsored affiliated packages will support *at least* the sunpy LTS version at the time of their release.

For dependencies only needed to run our tests we will support versions released in the prior 12 months to the current date.

.. _NEP-0029: https://numpy.org/neps/nep-0029-deprecation_policy.html

Extra Packages and Groups
=========================

sunpy has a number of optional dependencies that are not required to run the core functionality of the package.
These dependencies are grouped into "extras" that can be installed with the package.

We have to three categories of extras:

* subpackages : These extra groups are named after each subpackage and list the requirements that enable the core functionality of that subpackage.
* optional : These extra groups are named after package or feature that they enable. These are not needed to import sunpy or any of its subpackages, but are needed to use some of the functionality.
  It is expected that subpackages will handle their import manually with a nice import error message if the package is not installed.
* groupings : These are the groups that are used to install multiple optional dependencies at once.
  These range from "all" which is every optional dependency to "core" which only the dependencies that are needed to run the core functionality of sunpy.
  Also includes the "test", "docs" and "dev" groups.
