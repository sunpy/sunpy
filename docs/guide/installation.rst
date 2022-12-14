.. _guide_installing:

************
Installation
************

The SunPy project `maintains a range of packages <https://sunpy.org/project/affiliated>`__ that leverage the wider ecosystem of scientific Python packages for solar physics.

conda
=====
Full instructions for installing using conda are in :ref:`installing`.
This is the recommended way to install sunpy, because it creates a fresh Python environment that is independent from any other Python install or packages on your system.

Updating a package
------------------
You can update to the latest version of any package by running:

.. code-block:: bash

    $ conda update <package_name>

pip
===
This is for installing ``sunpy`` within a Python distribution, where ``pip`` has been used to install packages.
You will want to make sure you are using a `Python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__.

Once the environment active, to acquire a full ``sunpy`` installation:

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

.. note::

    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.
    Do **not** install ``sunpy`` or other Python packages using ``sudo``.
    This error implies you have an incorrectly configured virtual environment or it is not activated.
