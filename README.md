# [SunPy](http://sunpy.org)
[![Downloads](https://img.shields.io/pypi/dm/sunpy.svg)](https://pypi.python.org/pypi/sunpy/) 
[![Latest Version](https://img.shields.io/pypi/v/sunpy.svg)](https://pypi.python.org/pypi/sunpy/) 
[![Build Status](https://secure.travis-ci.org/sunpy/sunpy.svg)](http://travis-ci.org/sunpy/sunpy)
[![Build status](https://ci.appveyor.com/api/projects/status/xow461iejsjvp9vl?svg=true)](https://ci.appveyor.com/project/sunpy/sunpy)[![Coverage Status](https://coveralls.io/repos/sunpy/sunpy/badge.svg?branch=master)](https://coveralls.io/r/sunpy/sunpy?branch=master) [![Code Health](https://landscape.io/github/sunpy/sunpy/master/landscape.svg)](https://landscape.io/github/sunpy/sunpy/master)
[![Code Issues](https://www.quantifiedcode.com/api/v1/project/9edd3e28230840038713e1c7dc3eb141/badge.svg)](https://www.quantifiedcode.com/app/project/9edd3e28230840038713e1c7dc3eb141)[![Research software impact](http://depsy.org/api/package/pypi/sunpy/badge.svg)](http://depsy.org/package/python/sunpy)

SunPy is an open-source Python library for solar physics data analysis. See [sunpy.org](http://sunpy.org) for more information about the project.

For some examples of using SunPy see our [gallery](http://docs.sunpy.org/en/stable/generated/gallery/index.html).


Installation
------------

The recommended way to install SunPy is
with [conda](https://www.continuum.io/downloads). To install SunPy once conda is
installed run the following two commands:

    $ conda config --append channels conda-forge
    $ conda install sunpy


If you want to develop SunPy you will need to install from git. The best way to
do this is to create a new conda environment and install the git version of
SunPy in it:

    $ conda config --append channels conda-forge
    $ conda create -n sunpy-dev python sunpy hypothesis pytest-mock
    $ source activate sunpy-dev
    $ conda remove sunpy
    $ git clone https://github.com/sunpy/sunpy.git sunpy-git
    $ cd sunpy-git
    $ pip install -e .

For detailed installation instructions, see
the
[installation guide](http://docs.sunpy.org/en/latest/guide/installation/index.html) in
the SunPy docs.

Usage
-----

Here is a quick example of plotting an AIA image:

```python
>>> import sunpy.map
>>> from sunpy.data.sample import AIA_171_IMAGE
>>> aia = sunpy.map.Map(AIA_171_IMAGE)
>>> aia.peek()
```

Getting Help
------------

For more information or to ask questions about SunPy, check out:

 * [SunPy Documentation](http://docs.sunpy.org/en/latest/)
 * [SunPy Mailing List](https://groups.google.com/forum/#!forum/sunpy)
 * [SunPy Matrix Channel](https://riot.im/app/#/room/#sunpy-general:matrix.org)
 
Contributing
------------

If you would like to get involved, start by joining the
[SunPy mailing list](https://groups.google.com/forum/#!forum/sunpy)
and check out the [Developer's Guide](http://docs.sunpy.org/en/latest/dev.html) section
of the SunPy docs. Stop by our IRC chat room named #sunpy on irc.freenode.net
if you have any questions. Help is always welcome so let us know what you like
to work on, or check out the [issues page](https://github.com/sunpy/sunpy/issues)
for a list of some known outstanding items.

### Code of Conduct

When you are interacting with the SunPy community you are asked to follow
our [Code of Conduct](https://github.com/sunpy/sunpy/wiki/Code-of-Conduct).

