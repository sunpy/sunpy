=====
Linux
=====

Overview
--------

.. warning::
    We highly recommend that users use the Anaconda python distribution.
    Instructions are available here :ref:`anaconda_install`.

Almost all versions of Linux ship with a recent enough version
of Python, so it is unlikely that you will need to install Python yourself.
If you do not have python you can find the source at the
`Python official site <https://www.python.org/downloads/source/>`_.
Install using `sudo apt-get install python2.7-dev` or a similar command,
depending on your linux distro.
Anything like `python-dev` provides you the required Python development headers.

Ubuntu
------
On Ubuntu, most of the pre-reqs are available in the Ubuntu software repos and
can be installed using :command:`apt-get`: ::

    sudo apt-get install python-qt4
    sudo apt-get install git-core
    sudo apt-get python-numpy
    sudo apt-get python-scipy
    sudo apt-get python-matplotlib
    sudo apt-get update

Now we shall install pip.

Pip
^^^
Most Python distributions ship with a tool called
`easy_install <http://pypi.python.org/pypi/setuptools>`_
which assists with installing Python packages.

Although `easy_install`_ is capable of installing most of
the dependencies needed for SunPy itself, a more powerful tool called
`pip <http://pypi.python.org/pypi/pip>`__ provides a more flexible installation
(including support for uninstalling, upgrading, and installing from remote
sources such as GitHub) and should be used instead.

Use `easy_install`_ to install `pip`: ::

 sudo easy_install pip

You are now ready to :doc:`install SunPy and its dependencies <index>`.
