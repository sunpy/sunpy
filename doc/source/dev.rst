=================
Developer's Guide
=================

Getting Started
---------------
This article describes the guidelines to be followed by developers working on
SunPy.

Version Control
---------------

Source-code for SunPy is managed using the `Bazaar Distributed Version Control 
System (VCS) <http://bazaar.canonical.com/en/'>`_. Code branches are hosted on 
`Launchpad.net <http://launchpad.net/sunpy>`_, a free project hosting  website 
for Open-Source software.

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

Coding Standards
----------------
All code that is part of the SunPy project should follow The Syle Guide for 
Python (`PEP 8 <http://www.python.org/dev/peps/pep-0008/>`_).

Documentation
-------------

Code should be documented following the guidelines in `PEP 8 
<http://www.python.org/dev/peps/pep-0008/>`_ and `PEP 257 (Docstring 
conventions) <http://www.python.org/dev/peps/pep-0257/>`_. `Sphinx 
<http://sphinx.pocoo.org/>`_ is used to automatically generate documentation.

Files
-----

Each file should begin with a docstring describing the overall purpose of the 
module. Below that meta-tags containing author, license, email and credits 
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

Functions
---------

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
        
Testing
-------
Unit tests should be written as often as possible using unittest. See the 
Unit Testing section of Dive into Python 3 for more information about unit
testing in Python.

