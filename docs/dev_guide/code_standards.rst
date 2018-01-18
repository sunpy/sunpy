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
-  Documentation strings must be present for all public classes/methods/functions, and must follow the form outlined in the “Documentation Guidelines” document [tbc].
   Additionally, examples or tutorials in the package documentation are strongly recommended.
-  Unit tests are required for all public methods and functions, and should adhere to the standards set in the “Testing Guidelines”
   document [tbc].
-  C extensions will be allowed only when they provide a significant performance enhancement over pure python.
   When C extensions are used, the python interface must meet interface guidelines, and the use of Cython is strongly recommended.
-  If an external C library is needed, the source code for the library should be bundled with the SunPy core.
   Additionally, the package must be compatible with using a system-installed library.
-  Packages can include data in [path tbd] as long as it is less than about 100 kb.
   If the data exceeds this size, it should be hosted outside the source code repository and downloaded using the approved
   SunPy mechanism [tbc].
-  All persistent configuration should be stored using the functions in
   sunpy.config.
-  General utilities necessary for but not specific to the package should be placed in the .utils module.
   These utilities will be moved to the sunpy.utils module when the package is integrated into the core package. If a utility is already present in sunpy.utils, the package should always use that utility instead of re-implementing it in .utils.
   If changes or enhancements are required of that utility than a separate pull request should be presented to the community
-  Packages implementing many classes/functions not relevant to the component requested will not be accepted - the package should only
   include the required functionality and relevant extensions.
-  The use of short cryptic variable names is highly discouraged!
-  All code should follow The Style Guide for Python (`PEP 8 <http://www.python.org/dev/peps/pep-0008/>`_)
-  All code should follow the `coding style and conventions proposed by Astropy <http://docs.astropy.org/en/stable/development/codeguide.html#coding-style-conventions>`_.

PEP8 Rules
----------

Additionally, all code that goes in the project should be checked using `PyLint <http://www.logilab.org/card/pylint_manual>`_.
PyLint is an open source tool which analyzes Python code and checks for compliance with PEP8, as well as common coding errors and other potentially confusing or erroneous code statements.
Checking all submitted code helps to ensure some baseline level of quality and consistency for the code, and also helps to prevent potential problems from slipping through the cracks into the production code.

If you followed the installation instructions for developers, pylint should already be installed on your system.
To run PyLint on a file, simply call pylint from the command-line, passing in the name of the file you wish to check: ::

    pylint file.py

By default PyLint will print lines with potential problems along with a summary report.
To disable the summary report you can add either `-rn` or `--reports=no` to the command: ::

    pylint -rn file.py

Further, a paver task has been created so that all of the SunPy code can be checked at once: ::

    paver pylint

The output from PyLint will look something like: ::

 C: 87: Line too long (635/80)
 C:135: Line too long (98/80)
 R: 22:plot_fits: Too many local variables (22/15)
 R: 80:aia_color_table: Too many statements (59/50)
 W: 14: Unused import cm
 W: 16: Unused import Circle

Each line includes a line number, the category of the warning message, and a short description of the issue encountered.

The categories include:

* [R]efactor for a "good practice" metric violation
* [C]onvention for coding standard violation
* [W]arning for stylistic problems, or minor programming issues
* [E]rror for important programming issues (i.e. most probably bug)
* [F]atal for errors which prevented further processing

PyLint checks a wide range of different things so the first time you run PyLint on a file you will likely get a large number of warnings. In some cases the warnings will help you to spot coding mistakes or areas that could be improved with refactoring.
In other cases, however, the warning message may not apply and what you have there is exactly as it should be.
In these cases it is possible to silence PyLint for that line.
PyLint warning messages can be disabled at three different levels: globally (using a .pylintrc file), file-wide, and for a single line.

Porting Code from IDL
---------------------

As much IDL code for solar data analysis already exists in `SolarSoft <http://www.lmsal.com/solarsoft/>`__, it is expected that
many functions will need to be translated from IDL to Python for inclusion in SunPy.
The following describes the SunPy standards for how this should be done.

Naming conventions
~~~~~~~~~~~~~~~~~~

Much of the code in SolarSoft has its roots in Fortran.
Fortran is restrictive programming languages in terms of variable names therefore many variable and function names are short and cryptic.
The SunPy convention is to use easily understandable and more verbose names.
Many Python shells provide automatic code completion so the user experience is not likely to be impacted by longer names.
When converting IDL code, the code should be converted to a more Pythonic style (e.g. more verbose, clearer) when applicable unless the function is very well-known and already has an intuitive name.

Documentation
~~~~~~~~~~~~~

Copying the explanation (being careful to make modifications if things have changed) from the original IDL documentation, and including a
reference to the original IDL function is required.
While it may be informative to have the entire IDL documentation preserved in the Python code, this is discouraged.
