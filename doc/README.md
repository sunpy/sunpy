All of the SunPy documentation is contained in the doc/source folder and code
comments. To generate the documentation you must have Sphinx (as well as
Numpydoc and astropy-helpers) installed on your computer. Enter
the root folder and run:

  python setup.py build_sphinx -l

This will generate fresh HTML documentation for SunPy. Add the -o option to open the
docs in your browser after a successful build

  python setup.py build_sphinx -o

For more options run:

  python setup.py build_sphinx --help

For more information on how to use Sphinx, consult the
[Sphinx documentation](http://sphinx-doc.org).
