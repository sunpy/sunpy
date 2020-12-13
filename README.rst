********
`SunPy`_
********

|Latest Version| |codecov| |matrix| |Research software impact| |DOI| |Powered by NumFOCUS|

.. |Latest Version| image:: https://img.shields.io/pypi/v/sunpy.svg
   :target: https://pypi.python.org/pypi/sunpy/
.. |matrix| image:: https://img.shields.io/matrix/sunpy:openastronomy.org.svg?colorB=%23FE7900&label=Chat&logo=matrix&server_fqdn=openastronomy.modular.im
   :target: https://openastronomy.element.io/#/room/#sunpy:openastronomy.org
.. |codecov| image:: https://codecov.io/gh/sunpy/sunpy/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/sunpy/sunpy
.. |Research software impact| image:: http://depsy.org/api/package/pypi/sunpy/badge.svg
   :target: http://depsy.org/package/python/sunpy
.. |DOI| image:: https://zenodo.org/badge/2165383.svg
   :target: https://zenodo.org/badge/latestdoi/2165383
.. |Powered by NumFOCUS| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
   :target: https://numfocus.org
.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/sunpy/sunpy/master?filepath=examples

SunPy is an open-source Python library for Solar Physics data analysis and visualization.
Our homepage `SunPy`_ has more information about the project.

For some examples of using SunPy see our `gallery`_, to see the latest changes in SunPy see our `Changelog`_.

.. _SunPy: https://sunpy.org
.. _gallery: https://docs.sunpy.org/en/stable/generated/gallery/index.html
.. _Changelog: https://docs.sunpy.org/en/stable/whatsnew/changelog.html

Installation
============

The recommended way to install SunPy is with `miniconda`_.
To install SunPy once conda is installed run the following two commands:

.. code:: bash

    $ conda config --append channels conda-forge
    $ conda install sunpy

For detailed installation instructions, see the `installation guide`_ in the SunPy docs.

.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _installation guide: https://docs.sunpy.org/en/stable/guide/installation/index.html

Developing
==========

If you want to develop SunPy you will need to install from GitHub.
For detailed installation instructions, see `Development installation`_ in the SunPy docs.

Usage
=====

Here is a quick example of plotting an AIA image:

.. code:: python

    >>> import sunpy.map
    >>> from sunpy.data.sample import AIA_171_IMAGE
    >>> aia = sunpy.map.Map(AIA_171_IMAGE)
    >>> aia.peek()

Getting Help
============

For more information or to ask questions about SunPy, check out:

-  `SunPy Documentation`_
-  `SunPy Element Channel`_
-  `SunPy Mailing List`_

.. _SunPy Documentation: https://docs.sunpy.org/en/stable/
.. _SunPy Element Channel: https://app.element.io/#/room/#sunpy:openastronomy.org
.. _SunPy Mailing List: https://groups.google.com/forum/#!forum/sunpy

Contributing
============

|Open Source Helpers|

If you would like to get involved, start by joining the `SunPy mailing list`_ and check out the `Developers Guide`_ section of the SunPy docs.
Stop by our chat room `#sunpy:openastronomy.org`_ if you have any questions.
Help is always welcome so let us know what you like to work on, or check out the `issues page`_ for the list of known outstanding items.

For more information on contributing to SunPy, please read our `Newcomers' guide`_.

.. |Open Source Helpers| image:: https://www.codetriage.com/sunpy/sunpy/badges/users.svg
   :target: https://www.codetriage.com/sunpy/sunpy

.. _SunPy mailing list: https://groups.google.com/forum/#!forum/sunpy
.. _Developers Guide: https://docs.sunpy.org/en/latest/dev_guide/index.html
.. _`#sunpy:openastronomy.org`: https://app.element.io/#/room/#sunpy:openastronomy.org
.. _issues page: https://github.com/sunpy/sunpy/issues
.. _Newcomers' guide: https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html
.. _Development installation:  https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html#setting-up-a-development-environment

Code of Conduct
===============

When you are interacting with the SunPy community you are asked to follow our `Code of Conduct`_.

.. _Code of Conduct: https://sunpy.org/coc
