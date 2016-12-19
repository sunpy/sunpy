====================
Testing new features
====================

New features are constantly being developed for SunPy by contributors from
around the world. These new features are normally submitted as a
:ref:`Pull Request <developer-s-guide>`. Due to the review process imposed
on pull requests to maintain the stability and consistency of the code base,
it may take some time for a new feature to get accepted into the `master`
branch, and a bit more to get into the next stable release version.

In order to make it easier for people to test out new features (and hopefully provide valuable
feedback to the authors) we've put together this tutorial.

This assumes you have installed sunpy using anaconda, though a similar result
can be obtained using `virtualenv`.

We start by creating a new environment, so your current one still works as
expected::

 conda create -n sunpy-featureX python=3.5 sunpy

This creates an environment called `sunpy-featureX` (change `featureX` for
something more meaningful in your case) and specifies python 3.5. You can change
this to an older version of Python if you need to. Then we can activate it as ::

 source activate sunpy-featureX

or in windows ::

  activate sunpy-featureX

We proceed to remove the installed version of sunpy (only in this environment)
::

 conda remove sunpy

and now we can proceed to install the development version by cloning the git
repository ::

 cd to/a/place/I/want
 git clone https://github.com/sunpy/sunpy.git
 cd sunpy

That will download the latest development (`master) version of sunpy. You can
stop here if you just want and now you can test the bleeding edge of SunPy!
However, let's continue and `install` someone's new feature awaiting for approval
as a pull request. First `find out which pull request you want to test
<https://github.com/sunpy/sunpy/pulls>`_. Each pull request has a number like
`#1234`, use that number to get the feature you want to test as ::

 git fetch origin pull/1234/head:featureX
 git checkout featureX

Remember, change the `1234` by the pull request number, and the `featurex` by
something meaningful to you.

Next step is to install sunpy in your new environment ::

 python setup.py install

Once you've done so, you have the new feature ready for testing. Just start
your python and play with it. Once you are done with testing go back to your
original environment by ::

 source deactivate

or for windows, only ::

 deactivate

Remember, each time you start a new shell you will be using your default
environment, so you have to use reactivate your new feature environment.

**Final Note**

Please, once you have tested it, if you find it works as you expected, go to the
particular pull request discussion page and give your comments (or at least
thank the person who implemented it - it's always nice to know that what you do
is useful for others).
