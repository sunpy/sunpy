=================
Developer's Guide
=================

.. _dev-reference-label:

Developer's Guide Overview
--------------------------
This article describes the guidelines to be followed by developers working on
SunPy. You if you are thinking of contributing to SunPy please read the following
carefully.

Version Control
---------------

Source-code for SunPy is managed using `Git <http://git-scm.com>`_,
a Distributed Version Control system. Code branches are hosted on
`GitHub.com <http://github.com/sunpy/sunpy>`_, a free project hosting  website
for Open-Source software.

Creating Your Own Repo
^^^^^^^^^^^^^^^^^^^^^^

**Overview**

Each person contributing to SunPy should create their own code repository on
GitHub by forking the master repository or repo. All development is then done on that
fork, using topic branches to isolate work on different features. New
contributors can then initiate pull requests to have their code incorporated
into the SunPy master repository. Regular contributors can become members of the
`SunPy team <https://github.com/sunpy>`_ on GitHub. Code will be reviewed by regular
contributors and comments will usually be provided before code is accepted.

**Getting Started**

Creating your own repo on GitHub is easy to do. If you followed the SunPy installation
instructions you should already have git installed. Go ahead and create a free account
on `create an account on GitHub <https://github.com/signup/free>`_. Github has some
`great resources to help <https://help.github.com/>`_. Here is a quick overview of the
process.

**Adding an SSH key to GitHub**

Next, you need to tell GitHub who you are. In order to push any code to GitHub
you need to create a public SSH key and associate it with your GitHub account.
For instructions on how this is done, see the article on GitHub on
`Setting up git <http://help.github.com/set-up-git-redirect>`_ under
"Set Up SSH Keys". You only need to do this once, although if you plan to
work from multiple computers you will need to go through the process for each
computer you wish to work on. Once you have created your account and
associated a public SSH key it, you are ready to go.

**Using HTTPS**

If you do not fancy using SSH you can access GitHub using HTTP/HTTPS.
A few things to note.
Using HTTP only allows cloning of public repositories, while HTTPS allows cloning of private repositories but also allows you to have push access.
This way you can type in your username and password to access your repositories.

**Identifying yourself**

Begin by identifying yourself to git (so all of your commits have this information) and logging in to GitHub: ::

 git config --global user.name "Firstname Lastname"
 git config --global user.email "your_email@youremail.com"

**Forking SunPy**

Each contributor to SunPy has their own copy of the SunPy master repo. When
working on the code, changes are made to this copied repo, and only when the
changes are completed, and have been verified to work, are they pull requested back
to the upstream repo. GitHub provides a simple mechanism to setup your own
personal repo by providing an option to `fork a repository
<http://help.github.com/fork-a-repo/>`_. When you create a fork of a GitHub
project, a copy of the repo will automatically be created for you, and a link
will be provided which you can use to download the code to your machine and
begin working on it.

To begin, fork the main SunPy repo on GitHub by clicking on the `Fork` button
on the `SunPy project page <https://github.com/sunpy/sunpy>`_

Next, you need to download the forked repository. Clone the fork to your
local machine, edit and run: ::

 git clone git@github.com:your_username/sunpy.git

or: ::

 git clone http://github.com/sunpy/sunpy.git

By default your fork of the repo on GitHub is identified by the name `origin`.
In order to keep the fork up to date with the main repo, it is useful to add it
as a `remote` in git: ::

 git remote add upstream https://github.com/sunpy/sunpy.git

To stay up to date you can grab the latest changes to the SunPy master using
the commands: ::

 git pull upstream master

This will merge the upstream code automatically with your code so you don't need to worry
about it overwriting your changes. After running either of these commands,
your local copy of your personal repo is just a copy of the main repo.
This is the same procedure that you will use in the future to keep yourself synchronised with the
main repo. To make sure everything is setup correctly, let's make some changes
to our personal local repo and push those to our personal repo on GitHub. Go ahead and modify one
of the files, or create a new file (and then run :command:`git add`).

Commit and push the changes to GitHub: ::

 git commit -a -m "My first commit"
 git push

You local repo is now synced with GitHub and ahead of the main repo as it contains
your personal contribution. Remember to commit after you've done a unit of work (i.e.
often). This will make it easier for you (in the future) and everyone else to understand
what you are doing. Also make sure to make your commit statements clear and understandable.

**Installing SunPy**

In order to use the version of SunPy located  in your personal repository.
You need to install it using the `setup.py` script located in the top-level folder.
The `setup.py` script has several flags: ::
`develop` : Installs SunPy and builds all external libraries.
`build` or `build_ext`:  (Re)Builds the external libraries.
`clean --all`: Cleans all build files

Use the `setup.py` script like so: ::

 sudo python setup.py develop

If you are interested in having different versions of sunpy in your
machine and you want to switch from one to another you could use
virtual environments. This is an easy task if you used conda as your
python package manager.

After a standard conda installation, :ref:`assuming you have also installed
the latest stable version of sunpy <main-install>`, you then proceed to create a new environment
as::

 conda create -n sunpy-dev python=2.7 sunpy

This will create a new environment called `sunpy-dev` with all of the
dependencies needed by sunpy. We the proceed to change to the new
environment::

 source activate sunpy-dev

Then we need to remove the stable version from this environment ::

 conda remove sunpy

to then install the version in your git repository ::

 cd to/sunpy/git/repository
 python setup.py develop

At this stage you can use the development version in which you are
working on.
If you want to go back to the stable installation you can just change
the environment by ::

 source deactivate

**Conclusion**

That's it! You now have your own personal SunPy repo to develop on. You could
hack away at it to your heart's content, pushing changes to your fork on GitHub to share
with others and to ensure that you have a backup online.

But what about when you want to start contributing back to the main SunPy
repo? That is the topic of the next section.

Branches
^^^^^^^^

Developers should create topic branches within their repos for most of their
main coding. Every repo starts with a single branch called `master`, which
seldom needs to be used. Instead, work on any particular feature, bug, or
portion of the code is done in its own separate branch. This way changes on
any particular issue are isolated from other unrelated changes. Users can even
work on several different branches simultaneously.

To create a new branch run: ::

 git branch branchname

To switch to the new branch: ::

 git checkout branchname

(or alternatively, :command:`git checkout -b branchname` will accomplish
the above).

Developers should create new branches for the features they are working on.
When they have finished making changes and the code has been tested and
verified to be working well, the code can be merged back into the SunPy
repo. This is usually done through something called a pull request.

Example Workflow
^^^^^^^^^^^^^^^^

**Before we get started**

Here is an example workflow for a SunPy developer on any given day. Before
beginning this tutorial, follow the above instructions to grab a copy of the
SunPy repo.

**Grabbing other people's changes**

The first thing you want to do before you start coding anything new is to pull
in the latest code that others have written since you last did any coding. To
do this, run :command:`git pull`: ::

    git pull upstream master

This will ensure that you don't edit a file that has changed since your last pull
which will lead to merge conflicts later on.

**Code away**

Assuming there are no merge conflicts (which shouldn't happen unless two people
are working on the same part of the same file), then you are ready to begin
coding. If there are conflicts check out our conflicts section.

**Push your changes to GitHub**

As you code away on your local repo, you will need to keep git aware of what you are doing
and also your remote copy up to date.

To add a file, create the file then run: ::

    git add <yourfilename>

If you delete a file run: ::

    git rm <yourfilename>

To move a file: ::

    git mv <source> <destination>

To check to see if git is happy run: ::

    git status

which will give you a report of what has happened so far. Once you are at a good stopping point you should
"commit" your changes. This will provide you an opportunity to describe what you have done so far. To do this type: ::

    git commit -a -m "description of your changes"

After doing this you are ready to push your changes to your repo online with the command: ::

    git push

The local and remote copies of your repo are now synced.

**Contributing to the main repo**

Once you have made your desired changes, and committed and pushed your personal
branch, you need to decide whether or not to merge those changes back into the
main SunPy repo. If the changes you made are finished and have been tested and proven
stable (see the testing section below), then they can be merged into SunPy.
For now, lets assume that
your changes are complete and they are ready to be added to the main SunPy repo.
All contributed code to SunPy must be submitted as a "pull request". To do this go to the github
website and to your repo (remember to select the branch) then click on the "Pull
Request" button (in the upper right hand corner next to the Fork button which you've
used before). All initial pull requests must be made to the master branch unless they are a fix for specific version.
This will submit your code to a review. You will likely
receive some constructive comments on your code. To address these you can simply work
on your code and push those changes to your local repo. Those changes will be reflected
in your pull request. Once a member of
the SunPy dev team approves your pull request then your code will be
merged into the main SunPy repo
and your code will be part of the main SunPy code. Congratulations!

And that's it! It may seem like a lot at first but once you go through the
motions a few times it becomes very quick.

**Conflict resolution**

It may so happen that when you try to sync with the main repo there is a conflict error.
This means that someone else has been working on the same section of code
that you have. In such cases, the merge
command will issue a conflict warning and will then expect you do the merge
yourself. You can type: ::

   git mergetool

to go through the conflicts. This command will likely open some merging tools
which are already available on your computer. For example, on Mac OS X, it will open
FileMerge (if you have XCode installed). You can check on your progress by typing: ::

   git status

Once you are done, you should then commit your changes, in this case
the resolution of the conflict with: ::

   git commit -m "Resolved conflict between my and online version of file.py"

You can then proceed to push this change up to your branch.

Coding Standards
----------------
All code that is part of the SunPy project should follow The Style Guide for
Python (`PEP 8 <http://www.python.org/dev/peps/pep-0008/>`_) and
the `coding style and conventions proposed by Astropy
<https://astropy.readthedocs.org/en/stable/development/codeguide.html#coding-style-conventions>`_.
Additionally, all
code that goes in the trunk should be checked using `PyLint
<http://www.logilab.org/card/pylint_manual>`_. PyLint is an open source tool
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

(To be finished...)

Global Settings
---------------
SunPy makes use of a settings file (:file:`sunpyrc`). This file contains a
number of global settings such as where files should be downloaded by default
or the default format for displaying times. When developing new functionality
check this file and make use of the default values if appropriate or, if needed,
define a new value. More information can be found in :doc:`guide/customization`.

Documentation
-------------

All code must be documented. Undocumented code will not be accepted into SunPy.
Documentation should follow the guidelines in `PEP 8
<http://www.python.org/dev/peps/pep-0008/>`_ and `PEP 257 (Docstring
conventions) <http://www.python.org/dev/peps/pep-0257/>`_. Documentation for
modules, classes, and functions should follow the `NumPy/SciPy documentation
style guide
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_. We provide
an example of good documentation below or you can just browse some of SunPy code
itself for examples. All of the SunPy documentation (like this page!) is built by Sphinx
and must therefore adhere to Sphinx guidelines.

Sphinx
^^^^^^

**Overview**

`Sphinx <http://sphinx.pocoo.org/>`_ is a tool for generating high-quality
documentation in various formats (HTML, pdf, etc) and is especially well-suited
for documenting Python projects. Sphinx works by parsing files written using a
`a Mediawiki-like syntax
<http://docutils.sourceforge.net/docs/user/rst/quickstart.html>`_ called
`reStructuredText <http://docutils.sourceforge.net/rst.html>`_. In addition
to parsing static files of reStructuredText, Sphinx can also be told to parse
code comments. In fact, in addition to what you are reading right now, the
`Python documentation <http://www.python.org/doc/>`_ was also created using
Sphinx.

**Usage**

All of the SunPy documentation is contained in the ``doc/source`` folder and code
comments. To generate the documentation you must have Sphinx
(as well as Numpydoc and astropy-helpers) installed on your computer.
Enter the ``doc/source`` folder and run: ::

    make html

This will generate HTML documentation for SunPy. To clean up and delete the
generated documentation run: ::

    make clean

For more information on how to use Sphinx, consult the `Sphinx documentation
<http://sphinx.pocoo.org/contents.html>`_.

The rest of this section will describe how to document the SunPy code in order
to guarantee that well-formatted documentation will be created.

**doctest**

The example codes in the Guide section of the docs are configured with the Sphinx
`doctest extension <http://sphinx-doc.org/ext/doctest.html>`_.
This will test the example code to make sure it runs correctly, it can be executed
using: ::

  sphinx-build -t doctest -b doctest ./ ../build

from inside the ``doc/source`` folder.

Use of quantities and units
"""""""""""""""""""""""""""

Much code perform calculations using physical quantities.  SunPy uses astropy's
`quantities and units <http://docs.astropy.org/en/stable/units/index.html>`__
implementation to store, express and convert physical quantities. New classes
and functions should adhere to SunPy's `quantity and unit usage guidelines
<https://github.com/sunpy/sunpy-SEP/blob/master/SEP-0003.md>`__.  This document
sets out SunPy's reasons and requirements for the usage of quantities and
units.  Briefly, SunPy's `policy <https://github.com/sunpy/sunpy-SEP/blob/master/SEP-0003.md>`__
is that *all user-facing function/object arguments which accept physical
quantities as input **MUST** accept astropy quantities, and **ONLY** astropy
quantities*.

Developers should consult the
`Astropy Quantities and Units page <http://docs.astropy.org/en/stable/units/index.html>`__
for the latest updates on using quantities and units.  The `astropy tutorial on quantities and units
<http://www.astropy.org/astropy-tutorials/Quantities.html>`__ also provides useful examples on their
capabilities.

Astropy provides the decorator `~astropy.units.quantity_input` that
checks the units of the input arguments to a function against the
expected units of the argument.  We recommend using this decorator to
perform function argument unit checks.  The decorator ensures that the
units of the input to the function are convertible to that specified
by the decorator, for example ::

    import astropy.units as u
    @u.quantity_input(myangle=u.arcsec)
    def myfunction(myangle):
        return myangle**2

This function only accepts arguments that are convertible to arcseconds.
Therefore, ::

    >>> myangle(20 * u.degree)
    <Quantity 400.0 deg2>

returns the expected answer but ::

    >>> myangle(20 * u.km)

raises an error.

The following is an example of a use-facing function that returns the area of a
square, in units that are the square of the input length unit::

    @u.quantity_input(side_length=u.m)
    def get_area_of_square(side_length):
        """
        Compute the area of a square.

        Parameters
        ----------
        side_length : `~astropy.units.quantity.Quantity`
            Side length of the square

        Returns
        -------
        area : `~astropy.units.quantity.Quantity`
            Area of the square.
        """

        return (side_length ** 2)

This more advanced example shows how a private function that does not accept
quantities can be wrapped by a function that does::

    @u.quantity_input(side_length=u.m)
    def some_function(length):
        """
        Does something useful.

        Parameters
        ----------
        length : `~astropy.units.quantity.Quantity`
            A length.

        Returns
        -------
        length : `~astropy.units.quantity.Quantity`
            Another length
        """

        # the following function either
        # a] does not accept Quantities
        # b] is slow if using Quantities
        result = _private_wrapper_function(length.convert('meters').value)

        # now convert back to a quantity
        result = Quantity(result_meters, units_of_the_private_wrapper_function)

        return result

In this example, the non-user facing function *_private_wrapper_function* requires a numerical input in units of
meters, and returns a numerical output.  The developer knows that the result of *_private_wrapper_function* is in the
units *units_of_the_private_wrapper_function*, and sets the result of *some_function* to return the answer in those
units.


Examples
^^^^^^^^

Modules
"""""""

Each module or package should begin with a docstring describing its overall
purpose and functioning. Below that meta-tags containing author, license, email
and credits information may also be listed.

Example: ::

    """This is an example module comment.

    An explanation of the purpose of the module would go here and will appear
    in the generated documentation
    """
    #
    # TODO
    #  Developer notes and todo items can be listed here and will not be
    #  included in the documentation.
    #
    __authors__ = ["Keith Hughitt", "Steven Christe", "Jack Ireland", "Alex Young"]
    __email__ = "keith.hughitt@nasa.gov"
    __license__ = "xxx"

For details about what sections can be included, see the section on `documenting
modules
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_ in the
NumPy/SciPy style guide.

Classes
"""""""

Class docstrings should include a clear and concise docstring explaining the
overall purpose of the class, required and optional input parameters, and the
return value. Additionally, notes, references and examples are encouraged.

Example (:class:`sunpy.map.Map`) ::

    """
    Map(data, header)

    A spatially-aware data array based on the SolarSoft Map object

    Parameters
    ----------
    data : numpy.ndarray, list
        A 2d list or ndarray containing the map data
    header : dict
        A dictionary of the original image header tags

    Attributes
    ----------
    header : dict
        A dictionary representation of the image header
    date : datetime
        Image observation time
    det : str
        Detector name
    inst : str
        Instrument name
    meas : str, int
        Measurement name. For AIA this is the wavelength of image
    obs : str
        Observatory name
    r_sun : float
        Radius of the sun
    name : str
        Nickname for the image type (e.g. "AIA 171")
    center : dict
        X and Y coordinate for the center of the sun in arcseconds
    scale: dict
        Image scale along the x and y axes in arcseconds/pixel

    Examples
    --------
    >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    >>> aia.T
    Map([[ 0.3125,  1.    , -1.1875, ..., -0.625 ,  0.5625,  0.5   ],
    [-0.0625,  0.1875,  0.375 , ...,  0.0625,  0.0625, -0.125 ],
    [-0.125 , -0.8125, -0.5   , ..., -0.3125,  0.5625,  0.4375],
    ...,
    [ 0.625 ,  0.625 , -0.125 , ...,  0.125 , -0.0625,  0.6875],
    [-0.625 , -0.625 , -0.625 , ...,  0.125 , -0.0625,  0.6875],
    [ 0.    ,  0.    , -1.1875, ...,  0.125 ,  0.    ,  0.6875]])
    >>> aia.header['cunit1']
    'arcsec'
    >>> aia.show()
    >>> import matplotlib.cm as cm
    >>> import matplotlib.colors as colors
    >>> aia.peek(cmap=cm.hot, norm=colors.Normalize(1, 2048))

    See Also
    --------
    numpy.ndarray Parent class for the Map object

    References
    ----------
    | http://docs.scipy.org/doc/numpy/reference/arrays.classes.html
    | http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    | http://www.scipy.org/Subclasses

    """

Functions
"""""""""

Functions should include a clear and concise docstring explaining the overall
purpose of the function, required and optional input parameters, and the return
value. Additionally, notes, references and examples are encouraged.

Example (`numpy.matlib.ones
<https://github.com/numpy/numpy/blob/master/numpy/matlib.py>`_): ::

    def ones(shape, dtype=None, order='C'):
        """
        Matrix of ones.

        Return a matrix of given shape and type, filled with ones.

        Parameters
        ----------
        shape : {sequence of ints, int}
            Shape of the matrix
        dtype : data-type, optional
            The desired data-type for the matrix, default is np.float64.
        order : {'C', 'F'}, optional
            Whether to store matrix in C- or Fortran-contiguous order,
            default is 'C'.

        Returns
        -------
        out : matrix
            Matrix of ones of given shape, dtype, and order.

        See Also
        --------
        ones : Array of ones.
        matlib.zeros : Zero matrix.

        Notes
        -----
        If `shape` has length one i.e. ``(N,)``, or is a scalar ``N``,
        `out` becomes a single row matrix of shape ``(1,N)``.

        Examples
        --------
        >>> np.matlib.ones((2,3))
        matrix([[ 1.,  1.,  1.],
                [ 1.,  1.,  1.]])

        >>> np.matlib.ones(2)
        matrix([[ 1.,  1.]])

        """
        a = ndarray.__new__(matrix, shape, dtype, order=order)
        a.fill(1)
        return a

For details about what sections can be included, see the section on `documenting
functions
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_ in the
NumPy/SciPy style guide.

Trouble-shooting
^^^^^^^^^^^^^^^^
Sphinx can be very particular about formatting, and the warnings and errors
outputted aren't always obvious.

Below are some commonly-encountered warning/error messages along with a
human-readable translation:

**WARNING: Duplicate explicit target name: "xxx".**

If you reference the same URL, etc more than once in the same document sphinx
will complain. To avoid, use double-underscores instead of single ones after
the URL.

**ERROR: Malformed table. Column span alignment problem at line offset n**

Make sure there is a space before and after each colon in your class and
function docs (e.g. attribute : type, instead of attribute: type). Also, for
some sections (e.g. Attributes) numpydoc seems to complain when a description
spans more than one line, particularly if it is the first attribute listed.

**WARNING: Block quote ends without a blank line; unexpected unindent.**

Lists should be indented one level from their parents.

**ERROR: Unkown target name: "xxx"**

In addition to legitimate errors of this type, this error will also occur when
variables have a trailing underscore, e.g., ``xxx_``.

**WARNING: Explicit markup ends without a blank line; unexpected unindent.**

This usually occurs when the text following a directive is wrapped to the next
line without properly indenting a multi-line text block.

**WARNING: toctree references unknown document '...'** /
**WARNING: toctree contains reference to nonexisting document**

This pair of errors is due to the way numpydoc scrapes class members.

Testing
-------

This is a brief tutorial on how to write and run SunPy unit tests. SunPy makes use
of the great package `pytest <http://pytest.org>` for all of its testing needs.

Writing a unit test
^^^^^^^^^^^^^^^^^^^

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
straightforward in pytest: use the decorator ``@pytest.mark.online`` to
mark a test function as needing an internet connection.

Writing a unit test for a figure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
Once you have confirmed that the only figure-test failures are anticipated ones,
remove the existing hash library (found at `sunpy/tests/figure_hashes.json`)
and then run the entire suite of SunPy tests.  Note that all figure tests will
fail since a new hash library needs to be built.  The test report will tell you
where the new hash library has been created, which you then copy to
`sunpy/tests/`.

Running unit tests
^^^^^^^^^^^^^^^^^^

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

.. Unit tests should be written as often as possible using `unittest
.. <http://docs.python.org/release/3.1.3/library/unittest.html>`_. See the
.. `Unit Testing section <http://diveintopython3.org/unit-testing.html>`_ of
.. Dive into Python 3 for more information about unit testing in Python.

.. SunPy uses `tox <http://tox.testrun.org/>`_ to automate testing with
.. multiple versions of Python. The test environments are isolated and thus
.. all dependencies will need to be built; this requires the build dependencies
.. of those Python packages to be present on the system. These call be installed
.. by calling `sudo aptitude build-dep python-numpy python-scipy python-matplotlib python-pyfits`
.. on a distribution that derives from Debian. `tox` itself it also required and
.. can be installed by `pip install tox` (pip is a part of `python-distribute`).

.. The tests can then be run by running `tox` in the project directory.
.. This will take a very long time on the first run because it will
.. have to build all dependencies. Subsequent runs will take significantly
.. less time.

When to write unit tests
^^^^^^^^^^^^^^^^^^^^^^^^
A rule of thumb for unit testing is to have at least one unit test per public
function.

Testing Your Code Before Committing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When you commit your changes and make a Pull Request to the main SunPy repo on
GitHub, your code will be tested by Travis CI to make sure that all the tests
pass and the documentation builds without any warnings. Before you commit your
code you should check that this is the case. There is a helper script in
`sunpy/tools/pre-commit.sh` that is designed to run these tests automatically
every time you run `git commit` to install it copy the file from
`sunpy/tools/pre-commit.sh` to `sunpy/.git/hooks/pre-commit`, you should also
check the script to make sure that it is configured properly for your system.

Continuous Integration
^^^^^^^^^^^^^^^^^^^^^^

SunPy makes use of the `Travis CI service <https://travis-ci.org/sunpy/sunpy>`_.
This service builds a version of SunPy and runs all the tests. It also integrates
with GitHub and will report the test results on any Pull Request when they are
submitted and when they are updated.

The Travis CI server not only builds SunPy from source, but currently it builds all
of SunPy's dependencies from source as well using pip, all of this behaviour is
specified in the .travis.yml file in the root of the SunPy repo.

New Functionality
"""""""""""""""""
For SunPy, we would encourage all developers to thoroughly `cover <http://en.wikipedia.org/wiki/Code_coverage>`_
their code by writing unit tests for each new function created.

Developers who want to take an aggressive approach to reducing bugs may even
wish to consider adopting a practice such as Test Drive Development (TDD)
whereby unit tests are written before any actual code is written. The tests
begin by failing, and then as they code is developed the user re-runs the
tests until all of them are passing.

Bugs discovered
"""""""""""""""
In addition to writing unit tests new functionality, it is also a good practice
to write a unit test each time a bug is found, and submit the unit test along
with the fix for the problem. This way we can ensure that the bug does not
re-emerge at a later time.
