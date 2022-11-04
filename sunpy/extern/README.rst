****************
``sunpy.extern``
****************

This sub-package contains third-party Python packages/modules that are required for parts of ``sunpy``.

This currently includes the following:

- `AppDirs`_ This provides the core config file handling for ``sunpy``'s configuration system.
- `distro`_ This provides the functionality of the deprecated `platform.linux_distribution()`, which is used for ``sysinfo``.
- `inflect`_ This provides a way to transform numbers into words, which is used for the ``attrs`` system.
  **THIS IS NEVER UPDATED TO AVOID A NEW DEPENDENCY**
- `parse`_ This provides the functionality required for ``Scraper`` to parse URLs which powers ``GenericClient``.

.. _AppDirs: https://github.com/ActiveState/appdirs
.. _distro: https://github.com/nir0s/distro
.. _inflect: https://github.com/jaraco/inflect
.. _parse: https://github.com/r1chardj0n3s/parse
