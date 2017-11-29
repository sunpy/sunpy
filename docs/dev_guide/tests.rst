Testing
=======

This is a brief tutorial on how to write and run SunPy unit tests. SunPy makes use
of the great package `pytest <http://pytest.org>` for all of its testing needs.

Writing a unit test
-------------------

Consider a simple module `stuff.py` that contains the simple function shown
below.::

   def double(x):
       return 2 * x

We can write a test case for this function by defining a new function
containing the test (or tests) we want to perform. Suppose we want to check
that the correct behaviour occurs when we pass a value of 5 to `double()`. We
would write the test function like this: ::

  def test_answer():
      assert double(5) == 10

There are two things to note here. Firstly, names of test cases should always
begin with `test_`. This is because `pytest` searches for test cases named this
way. Secondly, we use `assert` to assert our expectation of what the result of
the test should be. In this example, the test returns true and so the test
passes.

The example given above is one in which the function and test reside in the
same module. In SunPy, functions and tests are separated and the latter can be
found in the `tests` directory within the directory containing the module.
The convention is to have one test module per module, with the names for
the test modules being the same as those for the modules prefixed with
`test_`. For example, the modules `xml.py` and `multimethod.py` in `sunpy/util`
have corresponding test modules `test_xml.py` and `test_multimethod.py` in
`sunpy/util/tests`.

There are some tests for functions and methods in SunPy that require a
working connection to the internet. pytest is configured in a way that it
iterates over all tests that have been marked as *online* and checks if
there is an established connection to the internet. If there is none, the
test is skipped, otherwise it is run. Marking tests is pretty
straightforward in pytest: use the decorator ``@pytest.mark.remote_data`` to
mark a test function as needing an internet connection.

Writing a unit test for a figure
--------------------------------

You can write SunPy unit tests that test the generation of matplotlib figures
by adding the decorator `sunpy.tests.helpers.figure_test`.
Here is a simple example: ::

    import matplotlib.pyplot as plt
    from sunpy.tests.helpers import figure_test

    @figure_test
    def test_simple_plot():
        plt.plot([0,1])

The current figure at the end of the unit test, or an explicitly returned
figure, has its hash compared against an established hash library (more on
this below).  If the hashes do not match, the figure has changed, and thus
the test is considered to have failed.

All such tests are automatically marked with the pytest mark
`pytest.mark.figure`.  See the next section for how to use marks.

You will need to update the library of figure hashes after you create a new
figure test or after a figure has intentionally changed due to code improvement.
After you have confirmed that any conflicting hashes are associated with desired
changes in figures, copy the hash-library file listed at the end of the test
report to `sunpy/tests/`.  Be forewarned that the hash library will likely need
to be updated for multiple versions of Python.

Running unit tests
------------------

To find and run all the SunPy unit tests, simply run ::

  py.test

from the root of the SunPy tree (i.e. the directory containing `INSTALL.TXT`,
`sunpy`, `doc`, etc.). This will produce a lot of output and you'll probably
want to run only selected test modules at a time. This is done by specifying
the module on the command line, e.g.::

 py.test sunpy/util/tests/test_xml.py

for the tests for `sunpy.util.xml`.

To run only tests that been marked with a specific pytest mark using the
decorator ``@pytest.mark`` (see the section *Writing a unit test*), use the
following command (where ``MARK`` is the name of the mark)::

  py.test -k MARK

To exclude (i.e. skip all tests with a certain mark, use the following
code (where ``MARK`` is the name of the mark)::

  py.test -k-MARK

Note that pytest is configured to skip all tests with the mark *online* if
there is no connection to the internet. This cannot be circumvented, i.e.
it cannot be forced to run a test with the mark *online* if there is no
working internet connection (rename the mark to something else to call the test
function anyway).

To get more information about skipped and xfailed tests (xfail means a
test has passed although it has been marked as ``@pytest.mark.xfail``),
you need to use the option ``-rs`` for skipped tests and ``-rx`` for
xfailed tests, respectively. Or use ``-rxs`` for detailed information on
both skipped and xfailed tests.


When to write unit tests
------------------------

A rule of thumb for unit testing is to have at least one unit test per public
function.


Writing Doctests
----------------

Code examples in the documentation will also be run as tests, this helps to
validate that the documentation is accurate and upto date. SunPy uses the same
doctest system as astropy, so for information on writing doctests see
:ref:`astropy:doctests` in the astropy documentation.



Testing Your Code Before Committing
-----------------------------------

When you commit your changes and make a Pull Request to the main SunPy repo on
GitHub, your code will be tested by Travis CI to make sure that all the tests
pass and the documentation builds without any warnings. Before you commit your
code you should check that this is the case. There is a helper script in
`sunpy/tools/pre-commit.sh` that is designed to run these tests automatically
every time you run `git commit` to install it copy the file from
`sunpy/tools/pre-commit.sh` to `sunpy/.git/hooks/pre-commit`, you should also
check the script to make sure that it is configured properly for your system.

Continuous Integration
----------------------

SunPy makes use of the `Travis CI service <https://travis-ci.org/sunpy/sunpy>`_.
This service builds a version of SunPy and runs all the tests. It also integrates
with GitHub and will report the test results on any Pull Request when they are
submitted and when they are updated.

The Travis CI server not only builds SunPy from source, but currently it builds all
of SunPy's dependencies from source as well using pip, all of this behaviour is
specified in the .travis.yml file in the root of the SunPy repo.

New Functionality
-----------------

For SunPy, we would encourage all developers to thoroughly `cover <http://en.wikipedia.org/wiki/Code_coverage>`_
their code by writing unit tests for each new function created.

Developers who want to take an aggressive approach to reducing bugs may even
wish to consider adopting a practice such as Test Drive Development (TDD)
whereby unit tests are written before any actual code is written. The tests
begin by failing, and then as they code is developed the user re-runs the
tests until all of them are passing.

Bugs discovered
---------------

In addition to writing unit tests new functionality, it is also a good practice
to write a unit test each time a bug is found, and submit the unit test along
with the fix for the problem. This way we can ensure that the bug does not
re-emerge at a later time.
