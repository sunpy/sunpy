----------------
Testing tutorial
----------------

This is a brief tutorial on how to write and run SunPy unit tests. To work
with the contents of this document, you will need to install the 
`pytest <http://pytest.org>`_ package, which can be found on 
`PyPI <http://pypi.python.org/pypi>`_.

1. Writing a unit test
----------------------

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
found in the `sunpy/tests` directory, organised by module. The convention is
to have one test module per science module, with the names for the test modules
being the same as those for the science modules prefixed with `test_`. For 
example, the modules `util.py` and `multimethod.py` in `sunpy/util` have 
corresponding test modules `test_util.py` and `test_multimethod.py`.

2. Running unit tests
---------------------

To find and run all the SunPy unit tests, simply run ::

  py.test

from the root of the SunPy tree (i.e. the directory containing `INSTALL.TXT`,
`sunpy`, `doc`, etc.). This will produce a lot of output and you'll probably 
want to run only selected test modules at a time. This is done by specifying
the module on the command line, e.g.::

 py.test sunpy/tests/util/test_util.py

for the tests for `sunpy.util.util`.
