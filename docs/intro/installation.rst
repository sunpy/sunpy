.. _installing:

************
Installation
************

The SunPy project `maintains a range of packages <https://sunpy.org/project/affiliated>`__ that leverage the wider ecosystem of scientific Python packages for solar physics.

Installing Python
=================
These instructions will set you up with `miniforge <https://conda-forge.org/docs/user/introduction.html>`__, which makes it easy to install and manage Python packages.

To install the miniforge Python distribution follow `the miniforge installation instructions <https://github.com/conda-forge/miniforge#install>`__.

This makes installing ``sunpy`` and many other packages in the scientific Python ecosystem much easier and quicker.
It also provides many pre-compiled binaries that are not available on PyPI.

Installing a package
--------------------
To install ``sunpy``, launch a terminal (under a UNIX-like system) or miniforge Prompt (under Windows).

It is considered good practice to create a new virtual environment for each project you work on.
`Here are some instructions on setting up a conda virtual environment. <https://towardsdatascience.com/getting-started-with-python-environments-using-conda-32e9f2779307>`__

You will want to activate your virtual environment then run:

.. code-block:: bash

    $ conda install sunpy

This will install ``sunpy`` and all of its dependencies.

You can replace ``sunpy`` in the above install command with other packages such as``sunkit-instruments`` to install them.
For a list of some other Python packages for solar physics, `see the SunPy Affiliated Packages list <https://sunpy.org/project/affiliated>`__.

.. warning::

    Please do not mix ``pip`` and ``conda`` to install packages within your environments.
    There are no guarantees that this will work and it is likely to cause problems.

Updating a package
------------------
You can update to the latest version of any package by running:

.. code-block:: bash

    $ conda update <package_name>

Installing without conda
========================
This is for installing ``sunpy`` within a Python distribution, where ``pip`` has been used to install packages.
You will want to make sure you are using a `Python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__.

Once the environment active, to acquire a full ``sunpy`` installation:

.. code-block:: bash

    $ pip install "sunpy[all]"

.. note::

    For fast and consistent installs, we strive to provide pure Python wheels for all of our packages.
    The core ``sunpy`` package also provides *compiled* wheels for CPython on Linux (x86-64) and macOS (x86-64 and ARM64), containing an optional C extension (``sunpy.io.ana``).
    On these platforms ``pip`` will install the compiled wheel by default, while on other platforms it will install the pure Python wheel without ``sunpy.io.ana``.
    If you require this extension on other platforms, install with ``pip install --no-binary sunpy "sunpy[all]"`` and ``pip`` will try to compile the package from source.
    Installing from source requires a C compiler (e.g., ``gcc`` or ``clang``), which can typically be installed from your system package manager, if not already installed.
    This extension is not compatible with Windows.

If you have a reason to want a more minimal installation, you can install sunpy with no optional dependencies, however this means a lot of submodules will not import:

.. code-block:: bash

    $ pip install "sunpy"

It is possible to select which "extra" dependencies you want to install, if you know you only need certain submodules:

.. code-block:: bash

    $ pip install "sunpy[map,timeseries]"

The available options are: ``[asdf]``, ``[dask]``, ``[database]``, ``[image]``, ``[jpeg2000]``, ``[map]``, ``[net]``, ``[timeseries]``, ``[visualization]``.

.. note::

    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.
    Do **not** install ``sunpy`` or other Python packages using ``sudo``.
    This error implies you have an incorrectly configured virtual environment or it is not activated.

If you want to develop ``sunpy`` we would strongly recommend reading the `Newcomers' Guide <https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html>`__.
