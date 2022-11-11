.. _conda_for_dependencies:

********************************************
Using ``conda`` for Development Dependencies
********************************************

The ``conda`` package for ``sunpy`` specifies only those dependencies that are relevant for the typical user.
A power user or a developer will instead need the complete set of dependencies to use all optional features, build the documentation, and/or run the full test suite.
We provide a ``conda`` environment file (``sunpy-dev-env.yml``) in the base directory of the GitHub repository.
You can create a ``conda`` environment (here named ``sunpy-dev``) from this file:

.. code-block:: bash

    $ conda env create -n sunpy-dev -f sunpy-dev-env.yml
    $ conda activate sunpy-dev

Alternatively, if you want to add the complete set of dependencies to an existing ``conda`` environment, you can use this file to update that environment (here named ``existing-env``):

.. code-block:: bash

    $ conda env update -n existing-env -f sunpy-dev-env.yml

The above calls assume that you have already downloaded ``sunpy-dev-env.yml``, either by itself or along with the whole repository.
You can alternatively supply the URL for the file, as hosted on GitHub, to either of the above calls, e.g.:

.. code-block:: bash

    $ conda env create -n sunpy-dev -f https://raw.githubusercontent.com/sunpy/sunpy/main/sunpy-dev-env.yml

.. note::

    This ``conda`` environment file intentionally does not specify restrictions on release versions for any dependency.

This ``conda`` environment file specifies only the *dependencies* for ``sunpy``, and not ``sunpy`` itself.
Depending on your needs, continue with one of the following two ways to install ``sunpy`` in this ``conda`` environment.

Normal Installation of ``sunpy``
================================

If you do not plan to modify the code of ``sunpy`` itself, you can simply install ``sunpy`` via ``conda``:

.. code-block:: bash

    $ conda install sunpy

Since the ``conda`` environment already has the complete set of dependencies for ``sunpy``, this call should install only ``sunpy`` and no additional packages.

Editable Installation of ``sunpy``
==================================

If you plan to modify the code of ``sunpy`` itself, you will want to perform an "editable install" of your local repository, via ``pip``, so that Python will link to your local repository.
Normally it is discouraged to have an environment that mixes package installations via ``conda`` with package installations via ``pip`` because it can lead to environment states that confuse the ``conda`` solver.
That is the reason why our instructions for new developers recommends that dependencies be installed exclusively via ``pip``.
However, some of the dependencies in the complete set are difficult or even impossible to install via ``pip`` alone, yet are straightforward to install via ``conda``.

Using the above ``conda`` environment, combined with a little care, it is possible to minimize that chance for any ``conda``/``pip`` conflicts.
From the base directory of your local repository, install ``sunpy`` with some additional ``pip`` options:

.. code:: bash

    $ pip install --no-deps --no-build-isolation -e .

The ``--no-deps`` and ``--no-build-isolation`` options ensure that ``pip`` does not itself install any dependencies.
Since the ``conda`` environment is designed to already have the complete set of dependencies, the ``pip`` installation should succeed.

You now have a ``conda`` environment with an editable installation of ``sunpy`` and with (nearly) all dependencies managed by ``conda``.
As you install other packages in this environment, a package that depends on ``sunpy`` will trigger ``conda`` to install ``sunpy``.
**That is fine!**
This ``conda`` installation of ``sunpy`` will simply mask the ``pip`` installation of ``sunpy``.
All you need to do is to remove the ``conda`` installation with the ``--force`` option so that dependencies are left undisturbed:

.. code-block:: bash

    $ conda remove --force sunpy

Once the ``conda`` installation of ``sunpy`` is removed, the ``pip`` installation of ``sunpy`` will automatically be accessible again.

.. note::

    For those who use ``mamba`` instead of ``conda``, most ``conda`` commands can be translated by simply substituting "mamba" for "conda".  However, ``mamba remove`` does not support the ``--force`` option, so you do in fact have to call ``conda remove``.

As a tip, you can follow a similar procedure to incorporate editable installations of other packages (e.g., ``astropy``) in a ``conda`` environment.
You first install the package via ``conda`` to ensure its dependencies are present, then you remove the package alone without disturbing the dependencies, and finally you perform the editable install of the package from the base directory of your local repository:

.. code-block:: bash

    $ conda install astropy
    $ conda remove --force astropy
    $ pip install --no-deps --no-build-isolation -e .
