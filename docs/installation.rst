.. _installing:

************
Installation
************

The SunPy project `maintains a range of packages <https://sunpy.org/project/affiliated>`__ that leverage the wider ecosystem of scientific Python packages for solar physics.

Python distribution
===================
These instructions will set you up with `miniforge <https://conda-forge.org/docs/user/introduction.html>`__, which makes it easy to install and manage Python packages.

`To install the miniforge Python distribution follow these instructions <https://github.com/conda-forge/miniforge#install>`__.

Installing miniforge will set the default conda channel to ``conda-forge``, which is required for installing ``sunpy`` and many other packages in the scientific Python ecosystem.

Installing a package using miniforge
------------------------------------
To install ``sunpy``, launch a terminal (under a UNIX-like system) or miniforge Prompt (under Windows).

.. note::

    We strongly recommend using a `Python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__ or a `conda virtual environment. <https://towardsdatascience.com/getting-started-with-python-environments-using-conda-32e9f2779307>`__

You will want to activate your virtual environment then run:

.. code-block:: bash

    $ conda install sunpy

This will install ``sunpy`` and each of its dependencies.

You can replace ``sunpy`` with other packages such as ``ablog`` or ``sunkit-instruments``.
For the full range of packages available, `see the SunPy Affiliated Packages list <https://sunpy.org/project/affiliated>`__.

Updating a package
------------------
You can update to the latest version of any package by running:

.. code-block:: bash

    $ conda update <package_name>

Installing on top of an existing scientific Python environment
--------------------------------------------------------------
This section assumes you already have everything setup, whether that be conda or a Python virtual environment.
These commands are to be executed within these environments.

conda
^^^^^

If you want to install ``sunpy`` within a pre-existing conda environment, you will want to activate the virtual environment and run:

.. code-block:: bash

    $ conda install sunpy

This assumes that your pre-existing conda environment is already using the "conda-forge" channel.
If this is not the case, please install miniforge (using the instructions above).

pip
^^^
This is for installing ``sunpy`` within a scientific Python distribution or environment, where ``pip`` has been used to install packages.

To acquire a fully working ``sunpy`` installation:

.. code-block:: bash

    $ pip install "sunpy[all]"

.. note::

    We strive to provide binary wheels for all of our packages.
    If you are using a Python distribution or operating system that is missing a binary wheel.
    ``pip`` will try to compile the package from source and this is likely to fail without a C compiler (e.g., ``gcc`` or ``clang``).
    Getting the compiler either from your system package manager or XCode should address this.

If you have a reason to want a more minimal installation, you can install sunpy with no optional dependencies, however this means a lot of submodules will not import:

.. code-block:: bash

    $ pip install "sunpy"

It is possible to select which "extra" dependencies you want to install, if you know you only need certain submodules:

.. code-block:: bash

    $ pip install "sunpy[map,timeseries]"

The available options are: ``[asdf]``, ``[dask]``, ``[database]``, ``[image]``, ``[jpeg2000]``, ``[map]``, ``[net]``, ``[timeseries]``, ``[visualization]``.

If you want to develop ``sunpy`` we would strongly recommend reading the `Newcomers' Guide <https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html>`__.

.. warning::

    Please do not mix ``pip`` and ``conda`` to install packages within your environments.
    There are no guarantees that this will work and it is likely to cause problems.

.. note::
    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.
    Do **not** install ``sunpy`` or other Python packages using ``sudo``.
    This error implies you have an incorrectly configured virtual environment or it is not activated.
