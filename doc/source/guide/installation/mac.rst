========
Mac OS X
========

Python
------
To begin, download and install the `latest 2.x version of Python <http://python.org/download/>`_
using the Mac OS X installer available from `Python.org <http://python.org/>`_.

Homebrew
--------
Next, install and update `homebrew <http://mxcl.github.com/homebrew/>`_: ::

 /usr/bin/ruby -e "$(/usr/bin/curl -fsSL https://raw.github.com/mxcl/homebrew/master/Library/Contributions/install_homebrew.rb)"
 brew doctor

Using homebrew, install Qt and some of the other dependencies needed for compilation later on by pip: ::

 brew -v  install gfortran pkgconfig git openjpeg qt

Pip
---
Most Python distributions ship with a tool called `easy_install <http://pypi.python.org/pypi/setuptools>`_ 
which assists with installing Python packages.

Although `easy_install`_ is capable of installing most of
the dependencies needed for SunPy itself, a more powerful tool called `pip <http://pypi.python.org/pypi/pip>`_
provides a more flexible installation (including support for uninstalling, upgrading, and installing from
remote sources such as GitHub) and should be used instead.

To begin, use `easy_install`_ to install `pip`: ::

 sudo easy_install pip

Use to install the remaining SunPy dependencies: ::

 sudo pip install --upgrade distribute
 sudo pip install --upgrade ipython
 sudo pip install --upgrade numpy
 sudo pip install --upgrade scipy
 sudo pip install --upgrade pyfits
 sudo pip install --upgrade suds
 sudo pip install --upgrade pandas
 sudo pip install --upgrade matplotlib
 
Additionally, if you plan to work help with SunPy development, some additional dependencies are required: ::

 pip install pytest pylint paver tox sphinx numpydoc

All done!
