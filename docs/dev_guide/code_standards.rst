Coding Standards
================

All code that is part of the SunPy project should follow The Style Guide for
Python (`PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_) and
the `coding style and conventions proposed by Astropy
<http://docs.astropy.org/en/stable/development/codeguide.html#coding-style-conventions>`_.
Additionally, all
code that goes in the trunk should be checked using `PyLint
<https://www.logilab.org/card/pylint_manual>`_. PyLint is an open source tool
which analyzes Python code and checks for compliance with PEP8, as well as
common coding errors and other potentially confusing or erroneous code
statements. Checking the SunPy trunk code this helps to ensure some baseline
level of quality and consistency for the code, and also helps to prevent
potential problems from slipping through the cracks into the production code.

If you followed the installation instructions for devs, pylint should already be
installed on your system. To run PyLint on a file, simply call pylint from the
command-line, passing in the name of the file you wish to check: ::

    pylint file.py

By default PyLint will print lines with potential problems along
with a summary report. To disable the summary report you can add either `-rn`
or `--reports=no` to the command: ::

    pylint -rn file.py

Further, a paver task has been created so that all of the SunPy code can be
checked at once: ::

    paver pylint

The output from PyLint will look something like: ::

 C: 87: Line too long (635/80)
 C:135: Line too long (98/80)
 R: 22:plot_fits: Too many local variables (22/15)
 R: 80:aia_color_table: Too many statements (59/50)
 W: 14: Unused import cm
 W: 16: Unused import Circle

Each line includes a line number, the category of the warning message, and a
short description of the issue encountered.

The categories include:

* [R]efactor for a "good practice" metric violation
* [C]onvention for coding standard violation
* [W]arning for stylistic problems, or minor programming issues
* [E]rror for important programming issues (i.e. most probably bug)
* [F]atal for errors which prevented further processing

PyLint checks a wide range of different things so the first time you run PyLint
on a file you will likely get a large number of warnings. In some cases the
warnings will help you to spot coding mistakes or areas that could be improved
with refactoring. In other cases, however, the warning message may not apply
and what you have there is exactly as it should be. In these cases it is
possible to silence PyLint for that line. PyLint warning messages can be
disabled at three different levels: globally (using a .pylintrc file),
file-wide, and for a single line.
