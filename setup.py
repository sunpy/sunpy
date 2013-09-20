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
    from setuptools import find_packages
    #Crotate Module
    from distutils.core import Extension
    from os.path import dirname, join
    cwd = dirname(__file__)
    try:
        import numpy as np
    except ImportError:
        print("SunPy WARNING: NumPy must be installed first to build the C extension")

    if 'np' in locals():
        module = 'sunpy.image.Crotate'   # import this
        sourcefiles = [join(cwd, 'sunpy', 'image', 'src', 'rot_extn.c'),
                       join(cwd, 'sunpy', 'image', 'src', 'transform', 'aff_tr.c')]
        libs = ['m']
        # -ON for compile optimise
        gcc_args = ['-std=c99', '-O3']
        # gcc_args = ['-std=c99']

        # need *module* name here
        crotate = Extension(module,
                            sources = sourcefiles,
                            libraries = libs,
                            extra_compile_args = gcc_args,
                            include_dirs =
                            [np.get_include(), join(cwd, 'sunpy', 'image', 'src')]
                            )

    setup(
	author="Steven Christe, Matt Earnshaw,  Russell Hewett, Keith Hughitt, Jack Ireland, Florian Mayer, Stuart Mumford,  Albert Shih, David Perez-Suarez et. al",
        author_email="sunpy@googlegroups.com",
        classifiers=CLASSIFIERS,
        description=DOCLINES[0],
        download_url="http://www.sunpy.org/download/",
        # 2011/11/21: disabling for now to prevent paver warnings
        #extra_requires={
        #    "Plotman": ['PyQt4']
        #},
        install_requires=[
            'numpy>1.6.0',
            'astropy>=0.2.0',
            'scipy',
 #           'suds',
            'pandas>=0.10.0',
            'matplotlib>=1.1',
 #           'beautifulsoup4',
        ],
        license="BSD",
        long_description="\n".join(DOCLINES[2:]),
        maintainer="SunPy Developers",
        maintainer_email="sunpy@googlegroups.com",
        name="sunpy",
        packages=find_packages(),
        package_data={'': ['*.fits', '*.fit', 'sunpyrc']},
        platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
        provides=['sunpy'],
        url="http://www.sunpy.org/",
        use_2to3=True,
        version="0.3.1",
        ext_modules = [crotate] if 'crotate' in locals() else []
    )

if __name__ == '__main__':
    from distribute_setup import use_setuptools
    use_setuptools()
    from setuptools import setup
    install(setup)
