
SunPy
=====

SunPy is an open-source Python library for solar physics data analysis.

Installation
------------

To begin, install the following requirements:

 * [Python]([Python](http://www.python.org) (2.6+)
 * [PyFITS](http://www.stsci.edu/resources/software_hardware/pyfits)
 * [NumPy](http://numpy.scipy.org/)
 * [SciPy](http://www.scipy.org/)
 * [Matplotlib](http://matplotlib.sourceforge.net/) (1.0+)
 * [Suds](https://fedorahosted.org/suds)
 * [pandas](http://pandas.pydata.org/)

Next, use git to grab the latest version of SunPy:

    git clone git://github.com/sunpy/sunpy.git
    cd sunpy
    python setup.py install

Done!

For detailed installation instructions, see the [installation guide](http://www.sunpy.org/doc/guide/installation/index.html) 
in the SunPy docs.

Usage
-----

Here is a quick example of plotting an AIA image:

```python
>>> import sunpy
>>> import matplotlib.cm as cm
>>> m = sunpy.make_map(sunpy.AIA_171_IMAGE)
>>> m.show(cmap=cm.hot)
```

Getting Help
------------

For more information or to ask questions about SunPy, check out:

 * [SunPy Documentation](http://www.sunpy.org/doc/)
 * [SunPy Mailing List](https://groups.google.com/forum/#!forum/sunpy)
 * IRC: #sunpy on [freenode.net](http://webchat.freenode.net/)

Contributing
------------

If you would like to get involved, start by joining the 
[SunPy mailing list](https://groups.google.com/forum/#!forum/sunpy)
and check out the [Developer's Guide](http://www.sunpy.org/doc/dev.html) section 
of the SunPy docs. Stop by our IRC chat room named #sunpy on irc.freenode.net 
if you have any questions. Help is always welcome so let us know what you like 
to work on, or check out the [issues page](https://github.com/sunpy/sunpy/issues) 
for a list of some known outstanding items.


