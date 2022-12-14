.. _installing:

************
Installation
************

This is the first chapter in the sunpy tutorial, and by the end of it you should have a working installation of Python and sunpy.

Installing Python
=================
There are many ways to install Python, but even if you have Python installed somewhere on your computer we recommend following these instructions anyway.
That's because we will create a new Python environment.
As well as containing a Python installation, this environment provides an isolated place to install Python packages (like sunpy).

The package manager we'll be using is called ``conda``.
To install conda, follow the `the miniforge installation instructions <https://github.com/conda-forge/miniforge#install>`__.
This will install ``conda`` and automatically configure the default channel (a channel is a remote software repository) to be ``conda-forge``, which is where ``sunpy`` packages are available.

Installing sunpy
----------------
To install ``sunpy``, start by launching a terminal (under a UNIX-like system) or miniforge Prompt (under Windows).
Then create a and activate new virtual environment

.. code-block:: bash

    $ conda create --name sunpy
    $ conda activate sunpy

In this case the environment is named 'sunpy'.
Feel free to change this to a different environment name.
Now we have a fresh environment we can install sunpy:

.. code-block:: bash

    $ conda install sunpy

This will install ``sunpy`` and all of its dependencies.
If you want to install another package later, you can run ``conda install <package_name>``.

Now we've got a working installation of sunpy, in the next few chapters we'll look at some of the basic data structures sunpy uses for representing times, coordinates, and data with physical units.

Futher reading
--------------
For more info on different ways to install ``sunpy`` (beyond the recommended way described above), see :ref:`guide_installing`.
