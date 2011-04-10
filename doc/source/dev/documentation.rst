Documentation
-------------

Code should be documented following the guidelines in PEP 8 and PEP 257 
(Docstring conventions). Sphinx will be used to automatically generate 
documentation.

Files
-----

Each file should begin with a docstring describing the overall purpose of the module. Below that meta-tags containing author, license, email and credits information should be listed.

Example:

::
    """This is an example module comment.
     
    An explanation of the purpose of the module would go here and will appear in the generated documentation
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

Functions should include a clear and concise docstring explaining the overall purpose of the function, required and optional input parameters, and the return value. Additionally, notes, references and examples are encouraged.

Example (numpy.matlib.ones):

::

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

