.. _sunpy-topic-guide-installing:

*********************
Advanced Installation
*********************

.. warning::

    This page has advanced instructions for installing ``sunpy``.
    If you are new to Python or ``sunpy``, please read the :ref:`installation instructions <sunpy-tutorial-installing>` first.

The SunPy Project `maintains a range of affiliated packages <https://sunpy.org/project/affiliated>`__ that leverage the wider ecosystem of scientific Python packages for solar physics.
This page focuses on how to install the ``sunpy`` core package, but these instructions should apply to most of the affiliated packages as well.

conda
=====

Full instructions for installing using conda are in :ref:`sunpy-tutorial-installing`.
This is the recommended way to install sunpy, because it creates a fresh Python environment that is independent from any other Python install or packages on your system, and allows the installation of non-python packages.
The SunPy Project publishes many of it's packages to the `conda-forge <https://conda-forge.org/>`__ channel.
If you have an existing conda install, without the conda-forge channel added you can configure the conda-forge channel by following the `instructions in the conda-forge <https://conda-forge.org/docs/user/introduction.html#how-can-i-install-packages-from-conda-forge>`__ documentation.

Updating a conda package
------------------------

You can update to the latest version of any package by running:

.. code-block:: bash

    $ conda update <package_name>

pip
===

This is for installing sunpy within a Python environment, where ``pip`` has been used to install all previous packages.
You will want to make sure you are using a `Python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__.

Once the environment active, to acquire a full ``sunpy`` installation:

.. code-block:: bash

    $ pip install "sunpy[all]"


.. warning::

    If you get a ``PermissionError`` this means that you do not have the required administrative access to install new packages to your Python installation.
    Do **not** install ``sunpy`` or other Python packages using ``sudo``.
    This error implies you have an incorrectly configured virtual environment or it is not activated.

We strive to provide binary wheels for all of our packages.
If you are using a Python version or operating system that is missing a binary wheel,
``pip`` will try to compile the package from source and this is likely to fail without a C compiler (e.g., ``gcc`` or ``clang``).
Getting the compiler either from your system package manager or XCode (if you are using macOS) should address this.

If you have a reason to want a more minimal installation, you can install sunpy with no optional dependencies, however this means a lot of submodules will not import:

.. code-block:: bash

    $ pip install "sunpy"

It is possible to select which "extra" dependencies you want to install, if you know you only need certain submodules:

.. code-block:: bash

    $ pip install "sunpy[map,timeseries]"

The available options are: ``[asdf]``, ``[dask]``, ``[image]``, ``[jpeg2000]``, ``[map]``, ``[net]``, ``[timeseries]``, ``[visualization]``.

Updating a pip package
----------------------

You can update to the latest version of any package by running:

.. code-block:: bash

    $ pip install --upgrade <package_name>
