.. _testing:

******************
Testing Guidelines
******************

This section describes the testing framework and format standards for tests in SunPy.
Here we have heavily adapted the `Astropy version <https://docs.astropy.org/en/latest/development/testguide.html>`_, and **it is worth reading that link.**

Testing Framework
=================

The testing framework used by SunPy is the `pytest`_ framework, accessed through the ``pytest`` command.

.. _pytest: https://pytest.org/en/latest/

.. note::

    The ``pytest`` project was formerly called ``py.test``, and you may
    see the two spellings used interchangeably.

Testing Dependencies
---------------------

Since the testing dependencies are not actually required to install or use SunPy, they are not included in "install_requires" in "setup.cfg".

Developers who want to run the test suite will need to install the testing packages using pip::

    $ pip install -e .[tests]

If you want to see the current test dependencies, you check "extras_require" in "setup.cfg".

Running Tests
-------------

There are currently two different ways to invoke the SunPy tests.
However, we strongly suggest using ``tox`` as the default one.
Each method uses the widely-used ``pytest`` framework and are detailed below.

``tox``
^^^^^^^^

The primary method is to use `tox`_, which is a generic virtualenv management and test command line tool.
We have several environments within our "tox.ini" file and you can list them::

    $ tox -v -l

Then you can run any of them doing::

    $ tox -e <name of env>

This will create a test environment in ".tox" and build, install SunPy and runs the entire test suite.
This is the method that our continuous integration uses.

.. _tox: https://tox.readthedocs.io/en/latest/

``pytest``
^^^^^^^^^^

The test suite can be run directly from the native ``pytest`` command.
In this case, it is important for developers to be aware that they must manually rebuild any extensions by running ``python setup.py build_ext`` before testing.

To run the entire suite with ``pytest``::

    $ pytest

will use the settings in ``setup.cfg``.

If you want to run one specific test file::

    $ pytest sunpy/map/tests/test_mapbase.py

or one specific test in a test file::

    $ pytest sunpy/map/tests/test_mapbase.py::<test_name>

(This does not work with ``tox`` and is a known issue.)

If a test errors, you can use ``pdb`` to create a debugging session at the moment the test fails::

    $ pytest --pdb

Test coverage reports
---------------------

SunPy can use `pytest-cov`_  generate test coverage reports and settings are stored in ``setup.cfg``.
This plugin can be installed using `pip`_::

    $ pip install pytest-cov

To generate a test coverage report, use::

    $ pytest --cov ./sunpy

This will print to the terminal a report of line coverage of our test suite.
If you want to create a report in html, you can run::

    $ pytest --cov-report xml:cov.xml --cov ./sunpy
    $ coverage html

.. _pytest-cov: https://pypi.org/project/pytest-cov/

Running tests in parallel
-------------------------

It is possible to speed up SunPy's tests using the `pytest-xdist`_ plugin.
This plugin can be installed using `pip`_::

    pip install pytest-xdist

Once installed, tests can be run in parallel using the ``--parallel`` commandline option.
For example, to use 4 processes::

    $ tox -e <name of environment> -- -n=4

or::

    $ pytest -n 4 ./sunpy

.. _pytest-xdist: https://pypi.python.org/pypi/pytest-xdist
.. _pip: https://pypi.org/project/pip/

Writing tests
=============

``pytest`` has the following `test discovery rules <https://pytest.org/en/latest/goodpractices.html#conventions-for-python-test-discovery>`_::

 * ``test_*.py`` or ``*_test.py`` files
 * ``Test`` prefixed classes (without an ``__init__`` method)
 * ``test_`` prefixed functions and methods

We use the first one for our test files, ``test_*.py`` and we suggest that developers follow this.

A rule of thumb for unit testing is to have at least one unit test per public function.

Simple example
--------------

The following example shows a simple function and a test to test this
function::

    def func(x):
        """Add one to the argument."""
        return x + 1

    def test_answer():
        """Check the return value of func() for an example argument."""
        assert func(3) == 5

If we place this in a ``test.py`` file and then run::

    $ pytest test.py

The result is::

    ============================= test session starts ==============================
    python: platform darwin -- Python 3.8.3 -- pytest-3.2.0
    test object 1: /Users/username/tmp/test.py

    test.py F

    =================================== FAILURES ===================================
    _________________________________ test_answer __________________________________

        def test_answer():
    >       assert func(3) == 5
    E       assert 4 == 5
    E        +  where 4 = func(3)

    test.py:5: AssertionError
    =========================== 1 failed in 0.07 seconds ===========================

Sometimes the output from the test suite will have ``xfail`` meaning a test has passed although it has been marked as ``@pytest.mark.xfail``), or ``skipped`` meaing a test that has been skipped due to not meeting some condition (online and figure tests are the most common).

You need to use the option ``-rs`` for skipped tests and ``-rx`` for xfailed tests, respectively.
Or use ``-rxs`` for detailed information on both skipped and xfailed tests.

Where to put tests
------------------

Each package should include a suite of unit tests, covering as many of the public methods/functions as possible.
These tests should be included inside each package, e.g::

    sunpy/map/tests/

"tests" directories should contain an ``__init__.py`` file so that the tests can be imported.

Online Tests
------------

There are some tests for functions and methods in SunPy that require a working connection to the internet.
``pytest`` is configured in a way that it iterates over all tests that have been marked as ``pytest.mark.remote_data`` and checks if there is an established connection to the internet.
If there is none, the test is skipped, otherwise it is run.

Marking tests is pretty straightforward, use the decorator ``@pytest.mark.remote_data`` to mark a test function as needing an internet connection::

    @pytest.mark.remote_data
    def func(x):
        """Add one to the argument."""
        return x + 1

By default, no online tests are selected and so to run the online tests you have to::

    $ tox -e py37-online

or::

    $ pytest --remote-data=any

Tests that create files
-----------------------

Tests may often be run from directories where users do not have write permissions so tests which create files should always do so in temporary directories.
This can be done with the `pytest tmpdir function argument <https://pytest.org/en/latest/tmpdir.html>`_ or with Python's built-in `tempfile module
<https://docs.python.org/3/library/tempfile.html#module-tempfile>`_.

Tests that use test data
------------------------

We store test data in "sunpy/data/test" as long as it is less than about 100 kB.
These data should always be accessed via the :func:`sunpy.data.test.get_test_filepath` and :func:`sunpy.data.test.test_data_filenames` functions.
This way you can use them when you create a test.

You can also use our sample data but this will have to be marked as an online test (see above)::

    import sunpy.data.sample

    @pytest.mark.remote_data
    def func():
        """Returns the file path for the sample data."""
        return sunpy.data.sample.AIA_131_IMAGE

Generally we do not run the tests on our sample data, so only do this if you have a valid reason.

Figure unit tests
-----------------

You can write SunPy unit tests that test the generation of matplotlib figures by adding the decorator `sunpy.tests.helpers.figure_test`.
Here is a simple example::

    import matplotlib.pyplot as plt
    from sunpy.tests.helpers import figure_test

    @figure_test
    def test_simple_plot():
        plt.plot([0,1])

The current figure at the end of the unit test, or an explicitly returned figure, has its hash (currently ``SHA256``) compared against an established hash collection (more on this below).
If the hashes do not match, the figure has changed, and thus the test is considered to have failed.

To run the figure tests you need to be very careful, as any pixel that has changed, will change the hash.
In order to avoid changes due to different package versions, we recommend using tox::

    $ tox -e figure

This will ensure that any figures created are checked using the package versions that were used to create the original figure hashes.
Running this will create a folder, "figure_test_images", within your work folder ("<local clone location>/figure_test_images"), which is ignored by git.
Inside this folder will be all the images created, as well as a json file with the hashes of the figures created by the test run.
The current hashes are located within "sunpy/tests/figure_tests_env_py36.json" and this will be where you will need to update old hashes or create new figure entries if anything changes.

If you are adding a new figure test, you will also need to update this `repository <https://github.com/sunpy/sunpy-figure-tests>`__ that stores the current figure tests, which we use for a visual comparison of figure tests.

Writing doctests
----------------

Code examples in the documentation will also be run as tests and this helps to validate that the documentation is accurate and up to date.
SunPy uses the same system as Astropy, so for information on writing doctests see the astropy `documentation <https://docs.astropy.org/en/latest/development/testguide.html#writing-doctests>`_.

You do not have to do anything extra in order to run any documentation tests.
Within our ``setup.cfg`` file we have set default options for ``pytest``, such that you only need to run::

    $ pytest <file to test>

to run any documentation test.

Bugs discovered
---------------

In addition to writing unit tests new functionality, it is also a good practice to write a unit test each time a bug is found, and submit the unit test along with the fix for the problem.
This way we can ensure that the bug does not re-emerge at a later time.
