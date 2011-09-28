"""
SunPy: Python for Solar Physics

The SunPy project is an effort to create an open-source software library for
solar physics using the Python programming language.
"""

DOCLINES = __doc__.split("\n")

CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
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

def install(setup): #pylint: disable=W0621
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
        packages=['sunpy', 'sunpy.cm', 'sunpy.map',
                  'sunpy.map.sources', 'sunpy.sun', 'sunpy.util',
                  'sunpy.tests', 'sunpy.tests.map', 'sunpy.tests.net',
                  'sunpy.tests.util',
                  'sunpy.data', 'sunpy.data.sample', 'sunpy.net', 'sunpy.solwcs',
                  'sunpy.net.vso', 'sunpy.net.hek',
                  'sunpy.gui', 'sunpy.gui.ui', 'sunpy.gui.ui.mainwindow',
                  'sunpy.gui.ui.mainwindow.resources', 'sunpy.gui.ui.mainwindow.widgets'],
        platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
        #install_requires=[
        #    'NumPy',
        #    'PyFITS',
        #    'SciPy',
        #    'Matplotlib>=1.0',            
        #],
        #extra_requires={
        #    "Plotman": ['PyQt4'],
        #    "VSO/HEK": ['suds']
        #},
        url="http://www.sunpy.org/",
        version="0.1",
        package_data={'': ['*.fits']},
    )

if __name__ == '__main__':
    from distutils.core import setup
    install(setup)
