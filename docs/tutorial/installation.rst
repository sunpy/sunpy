.. _sunpy-tutorial-installing:

************
Installation
************

This is the first chapter in the sunpy tutorial, and by the end of it you should have a working installation of Python and ``sunpy``.
For further information and alternative methods for installing ``sunpy`` beyond the recommended approach outlined below, refer to :ref:`sunpy-topic-guide-installing`.

Installing Python
=================

There are many ways to install Python, but even if you have Python installed somewhere on your computer we recommend following these instructions anyway.
That's because we will create a new Python virtual environment.
As well as containing a Python installation, this virtual environment provides an isolated place to install Python packages (like ``sunpy``) without affecting any other current Python installation.
If you already have Python and ``conda`` working you can skip the next section.

`If you are using Anaconda, we recommend that you uninstall it as the default package channel(s) have a restrictive license which means you might not be able to use it for free <https://sunpy.org/posts/2024/2024-08-09-anaconda/>`__.
Instead, we recommend that you use miniforge which is a minimal installer that setups ``conda`` with the ``conda-forge`` channel, which is free to use for everyone.
If you are using miniforge, you can skip the next section

.. _sunpy-tutorial-installing-miniforge:

Installing miniforge (and conda)
================================

If you don't already have a Python installation then we recommend installing Python with `miniforge <https://github.com/conda-forge/miniforge/#miniforge>`__.
Miniforge will install ``conda`` and automatically configure the default channel (a channel is a remote software repository) to be ``conda-forge``, which is where ``sunpy`` is available.

You will want to follow the instructions for your operating system and architecture on the `miniforge install page <https://conda-forge.org/download/>`__.
This will involve downloading a script or executable file and running it.

.. note::

   sunpy packages are not yet available on conda-forge for aarch64 or ppc64le
   architectures on linux, you will have to install sunpy from pip on these
   platforms.

To check that it has installed correctly, open a terminal (or the miniforge Prompt on Windows) and run:

.. code-block:: bash

    $ conda --version
    $ conda list

Installing sunpy
================

To install ``sunpy``, start by launching a terminal (under a UNIX-like system) or the miniforge Prompt (under Windows).
Now we will create and activate a new virtual environment to install ``sunpy`` into:

.. code-block:: bash

    $ conda create --name sunpy
    $ conda activate sunpy

In this case the virtual environment is named 'sunpy'.
Feel free to change this to a different environment name.

The benefit of using a virtual environment is that it allows you to install packages without affecting any other Python installations or versions on your system.
This also means you can work on multiple projects (research or coding) with different package requirements without them interfering with each other.

.. dropdown:: Click here if you haven't installed miniforge
    :color: warning

    If you have installed miniforge or are using Anaconda you need to configure conda to get your packages from conda-forge as well as the defaults channel.

    SunPy no longer recommends using the defaults channel at all, see `this blog post <https://sunpy.org/posts/2024/2024-08-09-anaconda/>`__ for details as to why.
    Therefore, if you are using Anaconda or miniconda we would suggest you uninstall it and install miniforge in its place.

    We also appreciate this isn't going to be possible for everyone, so what follows is our best instructions for how to proceed if you are using miniconda or Anaconda.

    The commands you need to run to add conda-forge and make it the default location to install conda packages from are:

    .. code-block:: bash

        $ conda config --add channels conda-forge
        $ conda config --set channel_priority strict

    These commands are taken from the
    `conda-forge documentation <https://conda-forge.org/docs/user/introduction/#how-can-i-install-packages-from-conda-forge>`__.

    Running these commands affect all the environments in your conda installation, critically, including the base Anaconda environment.
    We highly recommend that you do not install new packages, upgrade packages or use your base environment.
    Instead create new environments for all your projects, as you are much less likely to run into any pitfalls while using `multiple channels <https://conda-forge.org/docs/user/tipsandtricks/#multiple-channels>`__ by doing this.

Now that we have a fresh virtual environment, we can proceed with installing ``sunpy``:

.. code-block:: bash

    $ conda install sunpy

This will install ``sunpy`` and all of its dependencies.

To ensure that ``sunpy`` was installed correctly, run the following command:

.. code-block:: bash

    $ conda list sunpy

This checks if ``sunpy`` was installed correctly.

If you want to install another package later, you can run ``conda install <package_name>``.

Now we've got a working installation of ``sunpy``, in the next few chapters we'll look at some of the basic data structures ``sunpy`` uses for representing times, coordinates, and data with physical units.
