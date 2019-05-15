"""
This module provides functions to retrieve system information.
"""
import datetime
import platform

__all__ = ['get_sys_dict', 'system_info']


def get_sys_dict():
    """
    Test which packages are installed on system.

    Returns
    -------
    `dict`
        A dictionary containing the programs and versions installed on this machine.
    """
    try:
        from sunpy.version import version as sunpy_version
    except ImportError:
        sunpy_version = 'Missing sunpy version.'

    # Dependencies
    try:
        from numpy import __version__ as numpy_version
    except ImportError:
        numpy_version = "NOT INSTALLED"

    try:
        from scipy import __version__ as scipy_version
    except ImportError:
        scipy_version = "NOT INSTALLED"

    try:
        from matplotlib import __version__ as matplotlib_version
    except ImportError:
        matplotlib_version = "NOT INSTALLED"

    try:
        from astropy import __version__ as astropy_version
    except ImportError:
        astropy_version = "NOT INSTALLED"

    try:
        from pandas import __version__ as pandas_version
    except ImportError:
        pandas_version = "NOT INSTALLED"

    try:
        from bs4 import __version__ as bs4_version
    except ImportError:
        bs4_version = "NOT INSTALLED"

    try:
        from PyQt4.QtCore import PYQT_VERSION_STR as pyqt4_version
    except ImportError:
        pyqt4_version = "NOT INSTALLED"

    try:
        from PyQt5.QtCore import PYQT_VERSION_STR as pyqt5_version
    except ImportError:
        pyqt5_version = "NOT INSTALLED"

    try:
        from zeep import __version__ as zeep_version
    except ImportError:
        zeep_version = "NOT INSTALLED"

    try:
        from sqlalchemy import __version__ as sqlalchemy_version
    except ImportError:
        sqlalchemy_version = "NOT INSTALLED"

    try:
        from parfive import __version__ as parfive_version
    except ImportError:
        parfive_version = "NOT INSTALLED"

    try:
        from drms import __version__ as drms_version
    except ImportError:
        drms_version = "NOT INSTALLED"

    sys_prop = {'Time': datetime.datetime.utcnow().strftime("%A, %d. %B %Y %I:%M%p UT"),
                'System': platform.system(), 'Processor': platform.processor(),
                'SunPy': sunpy_version,
                'Arch': platform.architecture()[0], "Python": platform.python_version(),
                'NumPy': numpy_version, 'PyQt5': pyqt5_version,
                'SciPy': scipy_version, 'matplotlib': matplotlib_version,
                'Astropy': astropy_version, 'Pandas': pandas_version,
                'beautifulsoup': bs4_version, 'PyQt4': pyqt4_version,
                'Zeep': zeep_version, 'Sqlalchemy': sqlalchemy_version,
                'parfive': parfive_version, 'drms': drms_version
                }
    return sys_prop


def system_info():
    """
    Takes dictionary from sys_info() and prints the contents in an attractive
    fashion.
    """
    sys_prop = get_sys_dict()

    # title
    print("==============================")
    print("SunPy Installation Information")
    print("==============================\n")

    # general properties
    print("#######")
    print("General")
    print("#######")
    # OS and architecture information

    for sys_info in ['Time', 'System', 'Processor', 'Arch', 'SunPy']:
        print('{0} : {1}'.format(sys_info, sys_prop[sys_info]))

    if sys_prop['System'] == "Linux":
        distro = " ".join(platform.linux_distribution())
        print("OS: {0} (Linux {1} {2})".format(distro, platform.release(), sys_prop['Processor']))
    elif sys_prop['System'] == "Darwin":
        print("OS: Mac OS X {0} ({1})".format(platform.mac_ver()[0], sys_prop['Processor']))
    elif sys_prop['System'] == "Windows":
        print("OS: Windows {0} {1} ({2})".format(platform.release(), platform.version(),
                                                 sys_prop['Processor']))
    else:
        print("Unknown OS ({0})".format(sys_prop['Processor']))

    print("\n")
    # required libraries
    print("##################")
    print("Required Libraries")
    print("##################")

    for sys_info in ['Python', 'NumPy', 'SciPy', 'matplotlib', 'Astropy', 'Pandas', 'parfive']:
        print('{0}: {1}'.format(sys_info, sys_prop[sys_info]))

    print("\n")
    # recommended
    print("#####################")
    print("Recommended Libraries")
    print("#####################")

    for sys_info in ['beautifulsoup', 'PyQt4', 'PyQt5', 'Zeep', 'Sqlalchemy', 'drms']:
        print('{0}: {1}'.format(sys_info, sys_prop[sys_info]))
