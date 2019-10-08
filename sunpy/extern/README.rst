************
sunpy.extern
************

This sub-package contains third-party Python packages/modules that are
required for some of SunPy's core functionality.

In particular, this currently includes for Python:

- `AppDirs`_ (Python 2 and 3 versions): This provides the core config file handling for SunPy's configuration system.
- `distro` _ This provides the functionality of the deprecated `platform.linux_distribution()`, which we need (and more)

To replace any of the other Python modules included in this package, simply remove them and update any imports in SunPy to import the system versions rather than the bundled copies.

.. _AppDirs: https://github.com/ActiveState/appdirs
.. _distro: https://github.com/nir0s/distro.git
