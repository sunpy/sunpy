.. _coding-standards:


Coding Standards
================

The purpose of the page is to describe the standards that are expected of all code in sunpy.
All potential developers should read and abide by the following standards.
Code which does not follow these standards closely will not be accepted.

Language Standard
-----------------

English is the default language for all documentation strings and inline commands.
Variables names should also be based on English words.
In addition, the standard for spelling is American English.

-  Packages must be compatible with Python 2.7, and 3.x (for 3.x compatibility, the six package should be used for new code).
-  The code *should* be importable with no dependencies other than the Python Standard Library, NumPy, SciPy, matplotlib, and
   components already required by SunPy.
   Adding dependencies to SunPy will be considered but are highly discouraged.
-  Documentation strings must be present for all public classes/methods/functions, and must follow the form outlined in the :ref:`Documentation page <docs_guidelines>`.
   Additionally, examples or tutorials in the package documentation are strongly recommended.
-  Unit tests are required for all public methods and functions, and should adhere to the standards set in the :ref:`Testing page <testing>`.
-  C extensions will be allowed only when they provide a significant performance enhancement over pure python.
   When C extensions are used, the python interface must meet interface guidelines, and the use of Cython is strongly recommended.
-  If an external C library is needed, the source code for the library should be bundled with the SunPy core.
   Additionally, the package must be compatible with using a system-installed library.
-  Any data used for tests can be included in `sunpy/data/test` as long as it is less than about 1 MB.
-  If there is larger data files that are required, that can be added to the `sample-data repository <https://github.com/sunpy/sample-data>`_.
   Ideally, this data should not be used for code test but can be within documentation.
-  All persistent configuration should be stored using the functions in
   sunpy.config.
-  General utilities necessary for but not specific to the package should be placed in the .utils module.
   These utilities will be moved to the sunpy.utils module when the package is integrated into the core package. If a utility is already present in sunpy.utils, the package should always use that utility instead of re-implementing it in .utils.
   If changes or enhancements are required of that utility than a separate pull request should be presented to the community
-  Packages implementing many classes/functions not relevant to the component requested will not be accepted - the package should only
   include the required functionality and relevant extensions.
-  The use of short cryptic variable names is highly discouraged!
-  All code should follow The Style Guide for Python (`PEP 8 <http://www.python.org/dev/peps/pep-0008/>`_).
-  All code should follow the `coding style and conventions proposed by Astropy <http://docs.astropy.org/en/stable/development/codeguide.html#coding-style-conventions>`_.

PEP8 Rules
----------

Additionally, all code that goes in the project should be checked using one of the many PEP8 linters that are available (pep8, pylint, flake8 as some examples).
This is to ensure that the SunPy code follows The Style Guide for Python (`PEP 8 <http://www.python.org/dev/peps/pep-0008/>`_).

To ensure this, each pull request is checked for PEP8 issues by a bot called PEP8SPEAKS.
This way, people are aware of any potential issues with their submitted code.
