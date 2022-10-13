*********
``sunpy``
*********

|Latest Version| |DOI| |matrix| |codecov| |Binder| |Powered by NumFOCUS|

.. |Latest Version| image:: https://img.shields.io/pypi/v/sunpy.svg
   :target: https://pypi.python.org/pypi/sunpy/
.. |DOI| image:: https://zenodo.org/badge/2165383.svg
   :target: https://zenodo.org/badge/latestdoi/2165383
.. |matrix| image:: https://img.shields.io/matrix/sunpy:openastronomy.org.svg?colorB=%23FE7900&label=Chat&logo=matrix&server_fqdn=openastronomy.modular.im
   :target: https://openastronomy.element.io/#/room/#sunpy:openastronomy.org
.. |codecov| image:: https://codecov.io/gh/sunpy/sunpy/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/sunpy/sunpy
.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/sunpy/sunpy/main?filepath=examples
.. |Powered by NumFOCUS| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
   :target: https://numfocus.org

``sunpy`` is the core library of the `SunPy Project <https://sunpy.org/>`__, a Python software package which is designed to provide the fundamental tools for accessing, loading and interacting with solar physics data in Python.

For some examples of using ``sunpy`` see our `gallery <https://docs.sunpy.org/en/stable/generated/gallery/index.html>`__.
To see the latest changes in ``sunpy`` see our `changelog <https://docs.sunpy.org/en/stable/whatsnew/changelog.html>`__.

Installation
============

The recommended way to install ``sunpy`` is with `miniforge <https://github.com/conda-forge/miniforge#miniforge3>`__.
To install ``sunpy`` once conda is installed run the following command:

.. code:: bash

    $ conda install sunpy

For detailed installation instructions, see the `installation guide <https://docs.sunpy.org/en/stable/guide/installation.html>`__ in the ``sunpy`` docs.

Usage
=====

Here is a quick example of plotting an AIA image:

.. code:: python

   >>> import sunpy.map
   >>> from sunpy.data.sample import AIA_171_IMAGE
   >>>
   >>> import matplotlib.pyplot as plt
   >>>
   >>> aia = sunpy.map.Map(AIA_171_IMAGE)
   >>>
   >>> aia.plot()
   >>>
   >>> plt.show()

Getting Help
============

For more information or to ask questions about ``sunpy`` or any other SunPy library, check out:

-  `sunpy documentation <https://docs.sunpy.org/en/stable/>`__
-  `SunPy Chat <https://openastronomy.element.io/#/room/#sunpy:openastronomy.org>`__
-  `SunPy mailing list <https://groups.google.com/forum/#!forum/sunpy>`__

Acknowledging or Citing ``sunpy``
=================================

If you use ``sunpy`` in your scientific work, we would appreciate citing it in your publications.
The continued growth and development of ``sunpy`` is dependent on the community being aware of ``sunpy``.

Please see https://sunpy.org/about#acknowledging-or-citing-sunpy on how to do this.

Contributing
============

If you would like to get involved, start by joining the `SunPy Chat <https://openastronomy.element.io/#/room/#sunpy:openastronomy.org>`__ and check out our `Newcomers' guide <https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html>`__.
This will walk you through getting setup for contributing.

Code of Conduct
===============

When you are interacting with the SunPy community you are asked to follow our `Code of Conduct <https://sunpy.org/coc>`__.
