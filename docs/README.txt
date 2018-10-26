Documentation
=============

All of the SunPy documentation is contained in the `docs` folder and code
docstrings/comments. To generate the documentation you must have the packages
specified in `requirements/docs.txt` installed on your computer. Enter the root
folder and run:

    python setup.py build_docs -l

This will generate fresh HTML documentation for SunPy, and will completely
clean previous builds (including automodapi-generated files) before building
new ones.

Add the -o option to open the docs in your browser after a successful build:

    python setup.py build_docs -o

For more options run:

    python setup.py build_docs --help

For more information on how to use Sphinx, consult the
[Sphinx documentation](http://sphinx-doc.org).
