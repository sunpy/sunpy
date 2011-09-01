=================
Developer's Guide
=================

Overview
--------
This article describes the guidelines to be followed by developers working on
SunPy.

Version Control
---------------

Source-code for SunPy is managed using `Git <http://git-scm.com>`_, 
a Distributed Version Control system. Code branches are hosted on 
`GitHub.com <http://github.com/sunpy/sunpy>`_, a free project hosting  website 
for Open-Source software.

Creating Your Own Repo 
^^^^^^^^^^^^^^^^^^^^^^^^

**Overview**

Each person contributing to SunPy should create their own code repository on
GitHub by forking the master repository.All development is then done on that 
fork, using topic branches to isolate work on different features. New 
contributers can then initiate pull requests to have their code incorporated 
into the SunPy master repo. Regular contributers can become members of the 
SunPy team <https://github.com/sunpy>`_ on GitHub and push code directly to 
the master repo.

**Getting Started**

Creating your own repo on GitHub is easy to do. If you haven't already done so, 
`install git <http://git-scm.com/download>`_ and `create an account on 
GitHub <https://github.com/signup/free>`_.

**Adding an SSH key to GitHub**

Next, you need to tell GitHub who you are. In order to push any code to GitHub 
you need to create a public SSH key and associate it with your GitHub account. 
For instructions on how this is done, see the article on GitHub on 
`Setting up git <http://help.github.com/set-up-git-redirect>`_ under 
"Set Up SSH Keys". You only need to do this once, although if you plan to 
work from multiple computers you will need to go through the process for each 
computer you wish to work on. Once you have created your account and 
associated a public SSH key it, you are ready to go.

**Identifying yourself**

Begin by identifying yourself to GitHub and logging in to
GitHub: :: 

 git config --global user.name "Firstname Lastname"
 git config --global user.email "your_email@youremail.com"
 
**Forking SunPy**

Each contributer to SunPy has their own copy of the SunPy master repo. When
working on the code, changes are made to this copied repo, and only when the
changes are completed, and have been verified to work, are they pushed back
to the master repo. GitHub provides a simple mechanism to setup your own
personal repo by providing an option to `fork a repository 
<http://help.github.com/fork-a-repo/>`_. When you create a fork of a GitHub
project, a copy of the repo will automatically be created for you, and a link
will be provided which you can use to download the code to your machine and
begin working on it.

To begin, fork the main SunPy repo on GitHub by clicking on the `Fork` button 
on the `SunPy project page <https://github.com/sunpy/sunpy>`_

Next, you need to download the forked repository. Then clone the fork to your 
local machine, edit and run: ::

 git clone git@github.com:your_username/sunpy.git 
 
By default your fork of the repo on GitHub is identified by the name `origin`.
In order to keep the fork up to date with the main repo, it is useful to add it
as a `remote` in git: ::

 git remote add upstream git://github.com/sunpy/sunpy.git

To stay up to date you can grab the latest changes to the SunPy master using
the commands: ::

 git fetch upstream
 git merge upstream/master

Right now our local copy of your personal repo is just a copy of the main repo.
This is the same command you will use in the future to keep yourself syncronized with the
main repo. To make sure everything is setup correctly, let's make some changes
to our personal local repo and push those to our personal repo on GitHub. Go ahead and modify one
of the files, or create a new file (and then run :command:`git add`). 

Commit and push the changes to GitHub: ::

 git commit -a -m "My first commit"
 git push

You repo is now synced with GitHub and ahead of the main repo as it contains your personal contribution.

**Conclusion**

That's it! You now have your own personal SunPy repo to develop on. You could
hack away at it to your heart's content, pushing changes to GitHub to share
with others and to ensure that you have a backup online.

But what happens when you want to start contributing back to the main SunPy 
repo?

That is the topic of the next section.

Branches
^^^^^^^^^^^^^

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
verified to be working well, the code can be merged back into the SunPy master 
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
do this, run :command:`git fetch`: ::

    git fetch 
    
If no changes were made since the last time you worked on SunPy then you don't
need to do anything else and can begin coding again. If other people have pushed
code since you last worked on SunPy then these changes will be fetched and you
will need to merge them: ::

    git merge upstream/branch_name
    
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

To move a file, copy the file and then run a git rm and then a git add. To check to see if git is happy
run: ::

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
main repo. If the changes you made are finished and have been tested and proven
stable, then they can be merged into the main repo. If you are not finished with 
making your changes or broke some important functionality, then you will
probably want to wait before merging those changes. For now, lets assume that
your changes are complete and they are ready to be added to the main repo. The most
polite way to do this is to initiate a "pull request". To do this go to the github
website and to your repo (remember to select the branch) then click on the "Pull
Request" button (in the upper right hand corner next to the Fork button which you've
used before). This will submit your code to a review by the members of the SunPy dev
team. Once a member of the SunPy dev team approves it then your code is now part of
the main SunPy code. Congratulations!

And that's it! It may seem like a lot at first but once you go through the
motions a few times it becomes very quick.

**Conflict resolution**
It may so happen that when you try to sync with the main repo there is a conflict error.
This means that someone else has been working on the same section of code 
that you have. In such cases, the merge 
command will issue a conflict warning and will then expect you do the merge 
yourself. You can type: ::

   git mergetools

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
Python (`PEP 8 <http://www.python.org/dev/peps/pep-0008/>`_). Additionally, all
code that goes in the trunk should be checked using `PyLint 
<http://www.logilab.org/card/pylint_manual>`_. PyLint is an open source tool 
which analyzes Python code and checks for compliance with PEP8, as well as 
common coding errors and other potentially confusing or erroneous code 
statements. Checking the SunPy trunk code this helps to ensure some baseline
level of quality and consistency for the code, and also helps to prevent 
potential problems from slipping through the cracks into the production code.

PyLint can be installed using `easy_install`: ::

    easy_install pylint

To run PyLint on a file, simply call pylint from the command-line, passing in
the name of the file you wish to check: ::

    pylint file.py
    
By default PyLint will print a line of lines with potential problems along
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

Documentation
-------------

Code should be documented following the guidelines in `PEP 8 
<http://www.python.org/dev/peps/pep-0008/>`_ and `PEP 257 (Docstring 
conventions) <http://www.python.org/dev/peps/pep-0257/>`_. Documentation for 
modules, classes, and functions should follow the `NumPy/SciPy documentation 
style guide 
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_

Sphinx
^^^^^^

**Overview**

`Sphinx <http://sphinx.pocoo.org/>`_ is tool for generating high-quality 
documentation in various formats (HTML, pdf, etc) and is especially well-suited
for documenting Python projects. Sphinx works by parsing files written using a 
`a Mediawiki-like syntax 
<http://docutils.sourceforge.net/docs/user/rst/quickstart.html>`_ called 
`reStructuredText <http://docutils.sourceforge.net/rst.html>`_. In addition 
to parsing static files of reStructuredText, Sphinx can also be told to parse
code comments. In fact, in addition to what you are reading right now, the
`Python documenation <http://www.python.org/doc/>`_ was also created using
Sphinx.

**Usage**

All of the SunPy documentation is contained in the ``doc/source`` folder and code
comments. To generate the documentation you must have Sphinx (as well as Numpydoc) installed
on your computer (`easy_install sphinx` and `easy_install numpydoc`). Enter the ``doc/source`` folder and
run: ::

    make html

This will generate HTML documentation for SunPy.

Additionally, there is a `paver <http://paver.github.com/paver/>`_ command that
can be used to accomplish the same thing: ::

    paver build_sphinx

For more information on how to use Sphinx, consult the `Sphinx documentation 
<http://sphinx.pocoo.org/contents.html>`_.

The rest of this section will describe how to document the SunPy code in order
to guarantee that well-formatted documentation will be created.

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

Example (:class:`sunpy.map.BaseMap`) ::

    """
    BaseMap(data, header)
    
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
    >>> aia = sunpy.Map(sunpy.AIA_171_IMAGE)
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
    >>> aia.plot()
    >>> import matplotlib.cm as cm
    >>> import matplotlib.colors as colors
    >>> aia.plot(cmap=cm.hot, norm=colors.Normalize(1, 2048))
    
    See Also:
    ---------
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
        
Testing
-------
Unit tests should be written as often as possible using `unittest 
<http://docs.python.org/release/3.1.3/library/unittest.html>`_. See the 
`Unit Testing section <http://diveintopython3.org/unit-testing.html>`_ of 
Dive into Python 3 for more information about unit testing in Python.

SunPy uses `tox <http://tox.testrun.org/>`_ to automate testing with
multiple versions of Python. The test environments are isolated and thus
all dependencies will need to be built; this requires the build dependencies
of those Python packages to be present on the system. These call be installed
by calling `sudo aptitude build-dep python-numpy python-scipy python-matplotlib python-pyfits`
on a distribution that derives from Debian. `tox` itself it also required and
can be installed by `pip install tox` (pip is a part of `python-distribute`).

The tests can then be run by running `tox` in the project directory.
This will take a very long time on the first run because it will
have to build all dependencies. Subsequent runs will take significantly
less time.


Virtualenv
----------
`virtualenv <http://www.virtualenv.org/>`_ allows multiple isolated Python
environments to live on the same system. The `--no-site-packages` option
completely isolates it from the system Python installation; without it
packages installed on the system Python may also be used in the virtualenv.
