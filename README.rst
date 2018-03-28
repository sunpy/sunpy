`SunPy`_
========

|Latest Version| |Build Status| |Build status| |codecov| |Research software impact| |DOI| |Powered by NumFOCUS|

SunPy is an open-source Python library for solar physics data analysis.
See `sunpy.org`_ for more information about the project.

For some examples of using SunPy see our `gallery`_.

Installation
------------

The recommended way to install SunPy is with `conda`_.
To install SunPy once conda is installed run the following two commands:

.. code:: bash

    $ conda config --append channels conda-forge
    $ conda install sunpy

If you want to develop SunPy you will need to install from git.
The best way to do this is to create a new conda environment and install the git version of SunPy in it:

.. code:: bash

    $ conda config --append channels conda-forge
    $ conda create -n sunpy-dev python sunpy hypothesis pytest-mock
    $ source activate sunpy-dev
    $ conda remove sunpy
    $ git clone https://github.com/sunpy/sunpy.git sunpy-git
    $ cd sunpy-git
    $ pip install -e .

For detailed installation instructions, see the `installation guide`_ in the SunPy docs.

Usage
-----

Here is a quick example of plotting an AIA image:

.. code:: python

    >>> import sunpy.map
    >>> from sunpy.data.sample import AIA_171_IMAGE
    >>> aia = sunpy.map.Map(AIA_171_IMAGE)
    >>> aia.peek()

Getting Help
------------

For more information or to ask questions about SunPy, check out:

-  `SunPy Documentation`_
-  `SunPy Matrix Channel`_
-  `SunPy Mailing List`_

Contributing
------------

|Open Source Helpers|

If you would like to get involved, start by joining the `SunPy mailing list`_ and check out the `Developer’s Guide`_ section of the SunPy docs.
Stop by our chat room `#sunpy:matrix.org`_ if you have any questions.
Help is always welcome so let us know what you like to work on, or check out the `issues page`_ for the list of known outstanding items.

For more information on contributing to SunPy, please read our `contributing guide`_ or the `Newcomers guide`_.

Code of Conduct
~~~~~~~~~~~~~~~

When you are interacting with the SunPy community you are asked to
follow our `Code of Conduct`_.

.. |Latest Version| image:: https://img.shields.io/pypi/v/sunpy.svg
   :target: https://pypi.python.org/pypi/sunpy/
.. |Build Status| image:: https://secure.travis-ci.org/sunpy/sunpy.svg
   :target: http://travis-ci.org/sunpy/sunpy
.. |Build status| image:: https://ci.appveyor.com/api/projects/status/xow461iejsjvp9vl?svg=true
   :target: https://ci.appveyor.com/project/sunpy/sunpy
.. |codecov| image:: https://codecov.io/gh/sunpy/sunpy/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/sunpy/sunpy
.. |Research software impact| image:: http://depsy.org/api/package/pypi/sunpy/badge.svg
   :target: http://depsy.org/package/python/sunpy
.. |DOI| image:: https://zenodo.org/badge/2165383.svg
   :target: https://zenodo.org/badge/latestdoi/2165383
.. |Powered by NumFOCUS| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
   :target: http://numfocus.org
.. |Open Source Helpers| image:: https://www.codetriage.com/sunpy/sunpy/badges/users.svg
   :target: https://www.codetriage.com/sunpy/sunpy

.. _SunPy: http://sunpy.org
.. _sunpy.org: http://sunpy.org
.. _gallery: http://docs.sunpy.org/en/stable/generated/gallery/index.html
.. _conda: https://www.continuum.io/downloads
.. _installation guide: http://docs.sunpy.org/en/latest/guide/installation/index.html
.. _SunPy Documentation: http://docs.sunpy.org/
.. _SunPy Mailing List: https://groups.google.com/forum/#!forum/sunpy
.. _SunPy Matrix Channel: https://riot.im/app/#/room/#sunpy:matrix.org
.. _SunPy mailing list: https://groups.google.com/forum/#!forum/sunpy
.. _Developer’s Guide: http://docs.sunpy.org/en/latest/dev_guide/index.html
.. _`#sunpy:matrix.org`: https://riot.im/app/#/room/#sunpy:matrix.org
.. _issues page: https://github.com/sunpy/sunpy/issues
.. _contributing guide: https://github.com/sunpy/sunpy/blob/master/CONTRIBUTING.rst
.. _Newcomers guide: http://docs.sunpy.org/en/stable/dev_guide/newcomers.html
.. _Code of Conduct: http://docs.sunpy.org/en/stable/coc.html
