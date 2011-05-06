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

Creating Your Own Branch
^^^^^^^^^^^^^^^^^^^^^^^^

**Overview**

Each person contributing to SunPy should have their own development branch. This
is where most of the actual coding is done. Only when you have tested your
changes and are confident they are working properly and ready to be used by
others do they get pushed to the trunk.

Creating your own branch on Launchpad is pretty straight-forward. If you haven't
already done so, install Bazaar and `create an account on Launchpad 
<https://help.launchpad.net/YourAccount/NewAccount>`_. Next, you need to tell
Bazaar who you are. This is where the most difficult step of the process comes
in: in order to upload any code to Launchpad you need to create a public SSH key
and associate it with your launchpad account. For instructions on how this is 
done, see the article on Launchpad on `creating an SSH key pair 
<https://help.launchpad.net/YourAccount/CreatingAnSSHKeyPair>`_. Fortunately,
you will only need to do this once, although if you plan to work from multiple
computers you will need to go through the process for each computer you wish to
work on. Once you have created your account and associated a public SSH key it,
you are ready to go.

**Identifying yourself**

Begin by identifying yourself to Bazaar and logging in to
Launchpad: :: 

 bzr whoami "John Doe <john.doe@gmail.com>"
 bzr launchpad-login user-name
 
Next, grab a copy of the SunPy trunk: ::

 bzr branch lp:sunpy sunpy-trunk
 
This is the main branch where people's code gets merged together once it is
ready. For the most part, however, you won't be making changes to the trunk,
except to run an occasional `bzr pull` to keep it up to date.

**Branching from the trunk**

Next, let's branch from the trunk to create your own branch: ::

 bzr branch sunpy-trunk sunpy-dev
 
We could have also done `bzr branch lp:sunpy sunpy-dev`, but since we already
have a copy of the entire trunk on our computer we can just branch from that:
this is just a `small <http://www.joelonsoftware.com/items/2010/03/17.html>`_ 
`taste 
<http://doc.bazaar.canonical.com/migration/en/why-switch-to-bazaar.html>`_ of 
the benefits distrubted version control systems offer!

**Registering your branch on Launchpad**

Next, open up Launchpad and pull up the `SunPy code page 
<https://code.launchpad.net/sunpy>`_. On the right-hand side there should be
an a "Register a branch" option: click it. Give your branch a name that
identifies it as belonging to you (e.g. "jdoe") and choose "hosted" as the
branch type. Leave everything else as it is and hit "register branch" to create
your branch.

**Making changes and pushing them to Launchpad**

Once you have registered your branch on your Launch you will be taken to its
page. At the top of the page you will see a bzr command specifying how to push
to your new branch which should look something like `bzr push --use-existing 
lp:~john-doe/sunpy/jdoe`. Enter the directory for the branch you just created
with bzr and run this command to push the code to Launchpad: ::

 cd sunpy-dev
 bzr push --use-existing lp:~john-doe/sunpy/jdoe
 
If all goes well the changes should show up on Launchpad in a matter of seconds.
Right now our personal branch is just a copy of the trunk, which is not too
exciting. To make sure everything is setup correctly, let's make some changes
to our personal branch and push those to Launchpad. Go ahead and modify one
of the files, or create a new file (and then run `bzr add`) to the branch to
differentiate it from the trunk.

Commit and push the changes to Launchpad: ::

 bzr commit -m "My first commit"
 bzr push

Bazaar will remember the location you pushed to so you don't need to specify
the it again. Refresh the branch overview page on Launchpad and you should see
your changes reflected there once more.

**Conclusion**

That's it! You now have your own personal SunPy branch to develop on. You could
hack away at it to your heart's content, pushing changes to Launchpad to share
with others and to ensure that you have a backup online. Further, by issuing a
`bzr pull` command on the sunpy-trunk branch and then merging those changes into
your personal branch you can also stay up to date with changes made by others.

But what happens when you want to start contributing back to the SunPy trunk?
That is the topic of the next section.

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
programmer can use :command:`bzr pull` on the trunk followed by 
:command:`bzr merge` on their development branch to include any recent changes
made by other developers.

Example Workflow
^^^^^^^^^^^^^^^^

**Before we get started**

Here is an example workflow for a SunPy developer on any given day. Before
beginning this tutorial, follow the above instructions to grab a copy of the
SunPy trunk code and to set up your own branch. This tutorial assumes that you
have copies of both reposities on your computer. The personal branch will be
refered to as ``sunpy-dev`` and the trunk will be refered to as 
``sunpy-trunk``.

**Grabbing other people's changes**

The first thing you want to do before you start coding anything new is to pull
in the latest code that others have written since you last did any coding. To
do this, change directories to ``sunpy-trunk`` and run :command:`bzr pull`: ::

    bzr pull
    
If no changes were made since the last time you worked on SunPy then you don't
need to do anything else and can begin coding again. If other people have pushed
code since you last worked on SunPy then these changes will be fetched and you
will need to merge them into your development branch. To do this, enter your
personal branch and merge the changes from the trunk in: ::

    bzr merge ../sunpy-trunk
    
Make a quick commit in your branch to isolate the merge from any work you do.
    
**Code away**

Assuming there are no merge conflicts (which shouldn't happen unless two people
are working on the same part of the same file), then you are ready to begin
coding.

**Push your changes to Launchpad**

Once you have made your desired changes, and committed and pushed your personal
branch, you need to decide whether or not to merge those changes back into the
trunk. If the changes you made are finished and have been tested and proven
stable, then they can be merged into the trunk. If you are not finished making
with making your changes or broke some important functionality, then you will
probably want to wait before merging those changes into the trunk. For now, lets
assume that your changes are complete and they are ready to be added to the
trunk.

The first thing you will want to do is go into the trunk and run :command:`bzr 
pull` once more to see if any new changes have been made since you started 
coding: ::

    bzr pull

If there are new changes, then go ahead once more and merge those changes into
your personal branch and commit.

Next, change directories to the trunk and do a merge on your personal branch: ::

    bzr merge ../sunpy-dev
    
This will pull the changes you made into the trunk. Now all that remains is to
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

 TODO: Add PyLint instructions.

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
comments. To generate the documentation you must have Sphinx installed
on your computer (`easy_install sphinx`). Enter the ``doc/source`` folder and
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
    >>> map = sunpy.Map('doc/sample-data/AIA20110319_105400_0171.fits')
    >>> map.T
    Map([[ 0.3125,  1.    , -1.1875, ..., -0.625 ,  0.5625,  0.5   ],
    [-0.0625,  0.1875,  0.375 , ...,  0.0625,  0.0625, -0.125 ],
    [-0.125 , -0.8125, -0.5   , ..., -0.3125,  0.5625,  0.4375],
    ..., 
    [ 0.625 ,  0.625 , -0.125 , ...,  0.125 , -0.0625,  0.6875],
    [-0.625 , -0.625 , -0.625 , ...,  0.125 , -0.0625,  0.6875],
    [ 0.    ,  0.    , -1.1875, ...,  0.125 ,  0.    ,  0.6875]])
    >>> map.header['cunit1']
    'arcsec'
    >>> map.plot()
    >>> import matplotlib.cm as cm
    >>> import matplotlib.colors as colors
    >>> map.plot(cmap=cm.hot, norm=colors.Normalize(1, 2048))
    
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

