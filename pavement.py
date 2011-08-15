"""
SunPy: Python for Solar Physics

The SunPy project is an effort to create an open-source software library for
solar physics using the Python programming language.
"""
import os
import shutil
from paver.easy import *
import paver.doctools
from paver.setuputils import setup

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
DOCLINES = __doc__.split("\n")

CLASSIFIERS = [
    'Development Status :: 2 - Pre-Alpha',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Topic :: Software Development',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Physics',
    'Operating System :: Microsoft :: Windows',
    'Operating System :: POSIX',
    'Operating System :: Unix',
    'Operating System :: MacOS'
]
setup(
    author="Steven Christe, Keith Hughitt, Jack Ireland and Alex Young",
    author_email="keith.hughitt@nasa.gov",
    classifiers=CLASSIFIERS,
    description=DOCLINES[0],
    download_url="http://www.sunpy.org/download/",
    license="",
    long_description="\n".join(DOCLINES[2:]),
    maintainer="SunPy Developers",
    maintainer_email="sunpy@googlegroups.com",
    name="sunpy",
    packages=['sunpy', 'sunpy.cm', 'sunpy.dev', 'sunpy.map', 
              'sunpy.map.sources', 'sunpy.sun', 'sunpy.util'],
    platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
    url="http://www.sunpy.org/",
    version="0.01"
)

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
