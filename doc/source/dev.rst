=================
Developer's Guide
=================

Overview
--------
This article describes the guidelines to be followed by developers working on
SunPy.

Version Control
---------------

Source-code for SunPy is managed using the `Bazaar Distributed Version Control 
System (VCS) <http://bazaar.canonical.com/en/'>`_. Code branches are hosted on 
`Launchpad.net <http://launchpad.net/sunpy>`_, a free project hosting  website 
for Open-Source software.

Download the latest version of SunPy using Bazaar: ::

    bzr branch lp:sunpy

Collaboration
^^^^^^^^^^^^^

When multiple people are working on SunPy at the same time, the methods 
described by the `Team Collaboration/Distributed Development 
<http://doc.bazaar.canonical.com/latest/en/user-guide/distributed_intro.html>`_ 
article should be used as defined in the `Bazar User Guide 
<http://doc.bazaar.canonical.com/latest/en/user-guide/>`_.

Each developer should has his or her own personal development branch (e.g. 
"john-dev") where all of the main coding is done. When they have finished making
changes and the code has been tested and verified to be working well, the code 
can be merged back into the main trunk. When multiple developers are working on 
SunPy at the same time, care should be taken to avoid merging problems. Each 
user should keep a copy of the trunk on their own machine. Each day, the 
programmer can use "bzr pull" on the trunk followed by "bzr merge" on their 
development branch to include any recent changes made by other developers.

Example Workflow
^^^^^^^^^^^^^^^^

*Before we get started*
Here is an example workflow for a SunPy developer on any given day. Before
beginning this tutorial, follow the above instructions to grab a copy of the
SunPy trunk code and to set up your own branch. This tutorial assumes that you
have copies of both reposities on your computer. The personal branch will be
refered to simply as ``sunpy`` and the trunk will be refered to as 
``sunpy-trunk``.

*Grabbing other people's changes*
The first thing you want to do before you start coding anything new is to pull
in the latest code that others have written since you last did any coding. To
do this, change directories to ``sunpy-trunk`` and run bzr pull: ::

    bzr pull
    
If no changes were made since the last time you worked on SunPy then you don't
need to do anything else and can begin coding again. If other people have pushed
code since you last worked on SunPy then these changes will be fetched and you
will need to merge them into your development branch. To do this, enter your
personal branch and merge the changes from the trunk in: ::

    bzr merge ../sunpy-trunk
    
Make a quick commit in your branch to isolate the merge from any work you do.
    
*Code away*
Assuming there are no merge conflicts (which shouldn't happen unless two people
are working on the same part of the same file), then you are ready to begin
coding.

*Push your changes to Launchpad*
Once you have made your desired changes, and commited and pushed your personal
branch, you need to decide whether or not to merge those changes back into the
trunk. If the changes you made are finished and have been tested and proven
stable, then they can be merged into the trunk. If you are not finished making
with making your changes or broke some important functionality, then you will
probably want to wait before merging those changes into the trunk. For now, lets
assume that your changes are complete and they are ready to be added to the
trunk.

The first thing you will want to do is go into the trunk and run bzr pull once
more to see if any new changes have been made since you started coding: ::

    bzr pull

If there are new changes, then go ahead once more and merge those changes into
your personal branch and commit.

Next, change directories to the trunk and do a merge on your personal branch: ::

    bzr merge ../sunpy
    
This will pull the chnages you made into the trunk. Now all that remains is to
commit and push your changes back to Launchpad. While still in ``sunpy-trunk``,
run: ::

    bzr commit -m "description of your changes"
    bzr push

And that's it! It may seem like a lot at first but once you go through the
motions a few times it becomes very quick.


Coding Standards
----------------
All code that is part of the SunPy project should follow The Syle Guide for 
Python (`PEP 8 <http://www.python.org/dev/peps/pep-0008/>`_).

Documentation
-------------

Code should be documented following the guidelines in `PEP 8 
<http://www.python.org/dev/peps/pep-0008/>`_ and `PEP 257 (Docstring 
conventions) <http://www.python.org/dev/peps/pep-0257/>`_. Documentation for 
modules, classes, and functions should follow the `NumPy/SciPy documentation 
style guide 
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_


`Sphinx <http://sphinx.pocoo.org/>`_ is used to automatically generate 
documentation.

Module
^^^^^^

Each module should begin with a docstring describing its overall purpose and
functioning. Below that meta-tags containing author, license, email and credits 
information should be listed.

Example: ::

    """This is an example module comment.
     
    An explanation of the purpose of the module would go here and will appear 
    in the generated documentation
    """
    #
    # TODO
    #  Developer notes and todo items can be listed here
    #
    __author__ = "Keith Hughitt, Steven Christe, Jack Ireland and Alex Young"
    __email__ = "keith.hughitt@nasa.gov"
    __license__ = "MPL 1.0"

For details about what sections can be included, see the section on `documenting
modules 
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_ in the
NumPy/SciPy style guide.

Functions
^^^^^^^^^

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

