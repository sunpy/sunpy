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
# Packaging
#
DOCLINES = __doc__.split("\n")

CLASSIFIERS = [
    'Development Status :: 2 - Pre-Alpha',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: GNU General Public License (GPL)',
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
    packages=['sunpy'],
    platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
    url="http://www.sunpy.org/",
    version="0.01"
)

@task
@needs('html', 'generate_setup', 'setuptools.command.sdist')
def sdist():
    """Overrides sdist to make sure that our setup.py is generated."""
    shutil.rmtree('docs')
    pass

#
# Documentation
#
options(
    sphinx=Bunch(docroot='doc/source', builddir="_build")
)

@task
@needs('paver.doctools.html')
def html(options):
    """Build SunPy documentation"""
    sourcedir = 'doc/source/_build/html'
    destdir = 'docs'
    if os.path.exists(destdir):
        shutil.rmtree(destdir)
    shutil.move(sourcedir, destdir)
#
# Cleanup
#
@task
def clean():
    """Cleans up build files"""
    print("Removing build and sunpy.egg-info folders")
    for dir_ in ['dist', 'sunpy.egg-info']:
        if os.path.exists(dir_):
            shutil.rmtree(dir_)

