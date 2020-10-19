.. _coding-standards:

****************
Coding Standards
****************

The purpose of the page is to describe the standards that are expected of all the code in SunPy.
All potential developers should read and abide by the following standards.
Code which does not follow these standards closely will not be accepted.

We try to closely follow the coding style and conventions proposed by `Astropy <https://docs.astropy.org/en/stable/development/codeguide.html#coding-style-conventions>`_.

Language Standard
=================

* All code must be compatible with Python 3.7 and later.
  Usage of ``six``, ``__future__``, and ``2to3`` is no longer acceptable.

* The new Python 3 formatting style should be used (i.e.
  ``"{0:s}".format("spam")`` instead of ``"%s" % "spam"``).

* The core package and affiliated packages should be importable with no dependencies other than components already in SunPy, the `Python Standard Library <https://docs.python.org/3/library/index.html>`_, and packages already required by SunPy.
  Adding dependencies to SunPy will be considered but are highly discouraged.
  Such optional dependencies should be recorded in the ``setup.cfg`` file in the ``extras_require`` entry.

Coding Style/Conventions
========================

* The code will follow the standard `PEP8 Style Guide for Python Code <https://www.python.org/dev/peps/pep-0008/>`_.
  In particular, this includes using only 4 spaces for indentation, and never tabs.

* **Follow the existing coding style** within a file and avoid making changes that are purely stylistic.
  Please try to maintain the style when adding or modifying code.

* Following PEP8's recommendation, absolute imports are to be used in general.
  We allow relative imports within a module to avoid circular import chains.

* The ``import numpy as np``, ``import matplotlib as mpl``, and ``import matplotlib.pyplot as plt`` naming conventions should be used wherever relevant.
  ``from packagename import *`` should never be used (except in ``__init__.py``)

* Classes should either use direct variable access, or Python’s property mechanism for setting object instance variables.

* Classes should use the builtin :func:`super` function when making calls to methods in their super-class(es) unless there are specific reasons not to.
  :func:`super` should be used consistently in all subclasses since it does not work otherwise.

* Multiple inheritance should be avoided in general without good reason.

* ``__init__.py`` files for modules should not contain any significant implementation code. ``__init__.py`` can contain docstrings and code for organizing the module layout.

* General utilities necessary for but not specific to the package should be placed in the ``sunpy.utils`` module.

Formatting
==========

We enforce a minimum level of code style with our continuous intergration (the name is ``sunpy.sunpy (python_codestyle [linux]``).
This runs a tool called `pre-commit <https://pre-commit.com/>`__.

The settings and tools we use for the pre-commit can be found in the file `.pre-commit-config.yaml` at the root of the sunpy git repository.
Some of the checks are:
* Checks (but doesn't fix) various PEP8 issues with flake8.
* Sort all imports in any Python files with isort.
* Remove any unused variables or imports with autoflake.

We suggest you use "tox" (which is used to run the sunpy test suite) to run these tools without having to setup anything within your own Python virtual environment::

    $ tox -e codestyle

What you will see is this output (heavily condensed):

.. code-block:: bash

    codestyle create: /home/<USER>/GitHub/sunpy/.tox/codestyle
    codestyle run-test: commands[0] | pre-commit install-hooks
    codestyle run-test: commands[1] | pre-commit run --verbose --all-files --show-diff-on-failure
    flake8...................................................................Passed
    - hook id: flake8
    - duration: 1.35s

    0

    Check for case conflicts.................................................Passed
    - hook id: check-case-conflict
    - duration: 0.08s
    Trim Trailing Whitespace.................................................Failed
    - hook id: trailing-whitespace
    - duration: 0.08s
    - exit code: 1
    - files were modified by this hook

    Fixing docs/dev_guide/code_standards.rst

    pre-commit hook(s) made changes.
    If you are seeing this message in CI, reproduce locally with: `pre-commit run --all-files`.
    To run `pre-commit` as part of git workflow, use `pre-commit install`.
    All changes made by hooks:
    diff --git a/docs/dev_guide/code_standards.rst b/docs/dev_guide/code_standards.rst
    index bed700d90..c6b5df977 100644
    --- a/docs/dev_guide/code_standards.rst
    +++ b/docs/dev_guide/code_standards.rst
    @@ -59,6 +59,8 @@ Instead of installing this, you can use "tox" (which is used to run the sunpy te

        $ tox -e codestyle

    +What you will see
    +
    If you want to setup the pre-commit locally, you can do the following::

        $ pip install pre-commit
    diff --git a/docs/dev_guide/documentation.rst b/docs/dev_guide/documentation.rst
    index 5cd914047..b1017f77a 100644
    --- a/docs/dev_guide/documentation.rst
    +++ b/docs/dev_guide/documentation.rst
    @@ -39,9 +39,9 @@ If there are multiple code elements with the same name (e.g. ``peek()`` is a met

    .. code-block:: rst

    -    `GenericMap.peek` or `CompositeMap.peek`
    +    `.GenericMap.peek` or `.CompositeMap.peek`

    -These will show up as `GenericMap.peek` or `CompositeMap.peek`.
    +These will show up as `.GenericMap.peek` or `.CompositeMap.peek`.
    To still show only the last segment you can add a tilde as prefix:

    ERROR: InvocationError for command /home/nabil/GitHub/sunpy/.tox/codestyle/bin/pre-commit run --verbose --all-files --show-diff-on-failure (exited with code 1)
    ___________________________________________________________________________________________ summary ___________________________________________________________________________________________
    ERROR:   codestyle: commands failed

This will inform you of what checks failed and why, and what changes (if any) the command has made to your code.

If you want to setup the pre-commit locally, you can do the following::

    $ pip install pre-commit

Now you can do::

    $ pre-commit run --all-files

which will run the tools on all files in the sunpy git repository.
The pre-commit tools can change some of the files, but in other cases it will report problems that require manual correction.
If the pre-commit tool changes any files, they will show up as new changes that will need to be committed.

Automate
--------

Instead of running the pre-commit command each time you can install the git hook::

    $ pre-commit install

which installs a command to `.git/hooks/pre-commit` which will run these tools at the time you do ``git commit`` and means you don't have to run the first command each time.
We only suggest doing the install step if you are comfortable with git and the pre-commit tool.

Documentation and Testing
=========================

* American English is the default language for all documentation strings and inline commands.
  Variables names should also be based on English words.

* Documentation strings must be present for all public classes/methods/functions, and must follow the form outlined in the :ref:`docs_guidelines` page.
  Additionally, examples or tutorials in the package documentation are strongly recommended.

* Write usage examples in the docstrings of all classes and functions whenever possible.
  These examples should be short and simple to reproduce–users should be able to copy them verbatim and run them.
  These examples should, whenever possible, be in the :ref:`doctest <doctests>` format and will be executed as part of the test suite.

* Unit tests should be provided for as many public methods and functions as possible, and should adhere to the standards set in the :ref:`testing` document.

Data and Configuration
======================

* We store test data in ``sunpy/data/test`` as long as it is less than about 100 kB.
  These data should always be accessed via the :func:`sunpy.data.test.get_test_filepath` and :func:`sunpy.data.test.test_data_filenames` functions.

* We store data used for examples in the `sample-data repository <https://github.com/sunpy/sample-data>`_.
  This data should not be used for unit tests but can be within our documentation.

* All persistent configuration should use the :ref:`config` mechanism.
  Such configuration items should be placed at the top of the module or package that makes use of them, and supply a description sufficient for users to understand what the setting
  changes.

Standard output, warnings, and errors
=====================================

The built-in ``print(...)`` function should only be used for output that is explicitly requested by the user, for example ``print_header(...)`` or ``list_catalogs(...)``.
Any other standard output, warnings, and errors should follow these rules:

* For errors/exceptions, one should always use ``raise`` with one of the built-in exception classes, or a custom exception class.
  The nondescript ``Exception`` class should be avoided as much as possible, in favor of more specific exceptions (`IOError`, `ValueError`, etc.).

* For warnings, one should always use ``warnings.warn(message, warning_class)``.
  These get redirected to ``log.warning()`` by default, but one can still use the standard warning-catching mechanism and custom warning classes.
  The warning class should be either class:`~sunpy.utils.exceptions.SunPyUserWarning` or inherit from it.

Including C Code
================

* C extensions are only allowed when they provide a significant performance enhancement over pure Python, or a robust C library already exists to provided the needed functionality.

* The use of `Cython`_ is strongly recommended for C extensions.

* If a C extension has a dependency on an external C library, the source code for the library should be bundled with SunPy, provided the license for the C library is compatible with the SunPy license.
  Additionally, the package must be compatible with using a system-installed library in place of the library included in SunPy.

* In cases where C extensions are needed but `Cython`_ cannot be used, the `PEP 7 Style Guide for C Code <https://www.python.org/dev/peps/pep-0007/>`_ is recommended.

* C extensions (`Cython`_ or otherwise) should provide the necessary information for building the extension.

.. _Cython: https://cython.org/
