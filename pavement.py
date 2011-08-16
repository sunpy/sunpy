"""
SunPy: Python for Solar Physics

The SunPy project is an effort to create an open-source software library for
solar physics using the Python programming language.
"""
import os
import imp
import shutil
from paver.easy import *
import paver.doctools
from paver.setuputils import setup

# This is not pretty, but necessary
install = imp.load_source(
    'setup', os.path.join(os.path.dirname(__file__), 'setup.py')).install

#
# Options
#
options(
    deploy = Bunch(
        htmldir = path('doc/source/_build/html'),
        host = 'sipwork.org',
        hostpath = 'www/sunpy/doc'
    ),

    sphinx = Bunch(docroot='doc/source', builddir="_build"),
    pylint = Bunch(quiet=False)
)

#
# Packaging
#

install(setup)

@task
@needs('prepare_docs', 'setuptools.command.sdist')
def sdist():
    """Generated HTML docs and builds a tarball."""
    shutil.rmtree('doc/html')
    pass

#
# Documentation
#
@task
@needs('paver.doctools.html')
def prepare_docs(options):
    """Prepares the SunPy HTML documentation for packaging"""
    sourcedir = 'doc/source/_build/html'
    destdir = 'doc/html'
    if os.path.exists(destdir):
        shutil.rmtree(destdir)
    shutil.move(sourcedir, destdir)
    
@task
@needs('paver.doctools.html')
@cmdopts([('username=', 'u', 'Username')])
def deploy(options):
    """Update the docs on sunpy.org"""
    if "username" not in options:
        options.username = raw_input("Username: ")
    sh("rsync -avz -e ssh %s/ %s@%s:%s/" % (options.htmldir,
        options.username, options.host, options.hostpath))

#
# PyLint
#
@task
@cmdopts([('quiet', 'q', 'Run PyLint in quiet mode')])
def pylint(options):
    """Checks the code using PyLint"""
    from pylint import lint
    
    arguments = ['--rcfile=tools/pylint/pylintrc']
    
    if options.quiet:
        arguments.extend(["-rn"])
        
    arguments.extend(["sunpy/",  "tests/"])
    lint.Run(arguments)
    
#
# Cleanup
#
@task
@needs('paver.doctools.doc_clean')
def clean():
    """Cleans up build files"""
    print("Removing build files")
    for dir_ in ['doc/html', 'dist', 'sunpy.egg-info']:
        if os.path.exists(dir_):
            shutil.rmtree(dir_)
    for file_ in ['MANIFEST']:
        if os.path.exists(file_):
            os.remove(file_)
