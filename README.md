# [SunPy](http://sunpy.org)
[![Downloads](https://pypip.in/d/sunpy/badge.png)](https://pypi.python.org/pypi/sunpy/) [![Latest Version](https://pypip.in/v/sunpy/badge.png)](https://pypi.python.org/pypi/sunpy/) [![Build Status](https://secure.travis-ci.org/sunpy/sunpy.png)] (http://travis-ci.org/sunpy/sunpy)[![Build status](https://ci.appveyor.com/api/projects/status/xow461iejsjvp9vl?svg=true)](https://ci.appveyor.com/project/sunpy/sunpy)[![Coverage Status](https://coveralls.io/repos/sunpy/sunpy/badge.png?branch=master)](https://coveralls.io/r/sunpy/sunpy?branch=master) [![Code Health](https://landscape.io/github/sunpy/sunpy/master/landscape.png)](https://landscape.io/github/sunpy/sunpy/master)

SunPy is an open-source Python library for solar physics data analysis.

Installation
------------

To begin, install the following requirements:

 * [Python](http://www.python.org) (2.7+)
 * [Astropy](http://astropy.org) (1.0.0)
 * [NumPy](http://numpy.scipy.org/)
 * [SciPy](http://www.scipy.org/)
 * [Matplotlib](http://matplotlib.sourceforge.net/) (1.1+)
 * [Suds](https://fedorahosted.org/suds)
 * [pandas](http://pandas.pydata.org/) (0.10.0+)
 * [beautifulsoup4](http://www.crummy.com/software/BeautifulSoup/)
 * [sqlalchemy](http://www.sqlalchemy.org/)

Next, use git to grab the latest version of SunPy:

    git clone https://github.com/sunpy/sunpy.git
    cd sunpy
    python setup.py install

Done!

For detailed installation instructions, see the [installation guide](http://sunpy.readthedocs.org/en/latest/guide/installation/index.html)
in the SunPy docs.

Usage
-----

Here is a quick example of plotting an AIA image:

```python
>>> import sunpy.map
>>> import sunpy.data.sample
>>> import matplotlib.cm as cm
>>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
>>> aia.peek(cmap=cm.hot)
```

Getting Help
------------

For more information or to ask questions about SunPy, check out:

 * [SunPy Documentation](http://sunpy.readthedocs.org/en/latest/)
 * [SunPy Mailing List](https://groups.google.com/forum/#!forum/sunpy)
 * IRC: #sunpy on [freenode.net](http://webchat.freenode.net/)

Contributing
------------

If you would like to get involved, start by joining the
[SunPy mailing list](https://groups.google.com/forum/#!forum/sunpy)
and check out the [Developer's Guide](http://sunpy.readthedocs.org/en/latest/dev.html) section
of the SunPy docs. Stop by our IRC chat room named #sunpy on irc.freenode.net
if you have any questions. Help is always welcome so let us know what you like
to work on, or check out the [issues page](https://github.com/sunpy/sunpy/issues)
for a list of some known outstanding items.


