====================
Testing new features
====================

All the time new features are created to SunPy by contributors from around the
world. They are normally done as :ref:`Pull Request <developer-s-guide>`. These
pull request may take time to get part of the `master` branch, and a bit more to
get into the stable release version.

If you are brave enough and want to test new features, and provide feedback to
the authors, then you can follow up the steps below.

This assumes you have installed sunpy using anaconda, though a similar result
can be obtained using `virtualenv`.

We start creating a new environment, so your current one still works as
expected::

 conda create -n sunpy-featureX python=3.5 sunpy

This creates an environment called `sunpy-featureX` (change `featureX` for
something more meaningful in your case). Then we can activate it as ::

 source activate sunpy-featureX

We proceed to remove the installed version of sunpy (only in this environment)
::

 conda remove sunpy

and now we can proceed to install the development version by cloning the git
repository into a place where you want ::

 cd to/a/place/I/want
 git clone https://github.com/sunpy/sunpy.git
 cd sunpy

That will download the latest development version of sunpy, and if it's that
what you want to try then skip the next step. However, if you want to test
something that it's still as a pull request and awaiting for approval you have
to do as follows. First `find out which pull request you want to test
<https://github.com/sunpy/sunpy/pulls>`_. Each pull request has a number like
`#1234`, use that number to get the feature you want to test as ::

 git fetch origin pull/1234/head:featureX
 git checkout featureX

Remember, change the `1234` by the pull request number, and the `featurex` by
something meaningful to you.

Next step is to install sunpy in your new environment ::

 python setup.py install

Once you've done so, you have the new feature for testing in your computer, all
you need to do to use it is to do as usual, start your python and play with it.
Once you are done with testing such feature, go back to your default environment
by ::

 source deactivate

Remember, each time you start a new shell will be using your default, so you
have only to use ::

 source activate sunpy-featureX

to go back to the feature you installed. This will be just available from the
terminal you are using, so you could have two versions of sunpy at the same
time.

**Final Note**

Please, once you have tested it, if you find it works as you expected, go to the
particular pull request discussion page and give your comments (or at least
thank the person who implemented it - it's always nice to know that what you do
is useful for others).
