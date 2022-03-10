.. _dependency_versions:

*************************
Dependency Support Policy
*************************

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
