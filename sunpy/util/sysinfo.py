from __future__ import absolute_import, division, print_function

import platform
import datetime


__all__ = ['get_sys_dict', 'system_info']



def get_sys_dict():
    """
    Test which packages are installed on system.

    Returns
    -------
    sys_prop : `dict`
        A dictionary containing the programs and versions installed on this
        machine

    """

    try:
        from sunpy.version import version as sunpy_version
        from sunpy.version import githash as sunpy_git_description
    except ImportError:
        sunpy_version = 'Missing version.py; re-run setup.py'
        sunpy_git_description = 'N/A'

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
        from PyQt4.QtCore import PYQT_VERSION_STR as pyqt_version
    except ImportError:
        pyqt_version = "NOT INSTALLED"

    try:
        from suds import __version__ as suds_version
    except ImportError:
        suds_version = "NOT INSTALLED"

    try:
        from sqlalchemy import __version__ as sqlalchemy_version
    except ImportError:
        sqlalchemy_version = "NOT INSTALLED"

    try:
        from requests import __version__ as requests_version
    except ImportError:
        requests_version = "NOT INSTALLED"



    sys_prop = {'Time':datetime.datetime.utcnow().strftime("%A, %d. %B %Y %I:%M%p UT"),
                'System':platform.system(), 'Processor':platform.processor(),
                'SunPy':sunpy_version, 'SunPy_git':sunpy_git_description,
                'Arch':platform.architecture()[0], "Python":platform.python_version(),
                'NumPy':numpy_version,
                'SciPy':scipy_version, 'matplotlib':matplotlib_version,
                'Astropy':astropy_version, 'Pandas':pandas_version,
                'beautifulsoup':bs4_version, 'PyQt':pyqt_version,
                'SUDS':suds_version, 'Sqlalchemy':sqlalchemy_version, 'Requests':requests_version
                }
    return sys_prop

def system_info():
    """
    Takes dictionary from sys_info() and prints the contents in an attractive fashion.

    """
    sys_prop = get_sys_dict()

# title
    print("==========================================================")
    print(" SunPy Installation Information\n")
    print("==========================================================\n")


# general properties
    print("###########")
    print(" General")
    print("###########")
    # OS and architecture information

    for sys_info in ['Time', 'System', 'Processor', 'Arch', 'SunPy', 'SunPy_git']:
        print('{0} : {1}'.format(sys_info, sys_prop[sys_info]))

    if sys_prop['System'] == "Linux":
        distro = " ".join(platform.linux_distribution())
        print("OS: {0} (Linux {1} {2})".format(distro, platform.release(), sys_prop['Processor']))
    elif sys_prop['System'] == "Darwin":
        print("OS: Mac OS X {0} ({1})".format(platform.mac_ver()[0], sys_prop['Processor']))
    elif sys_prop['System'] == "Windows":
        print("OS: Windows {0} {1} ({2})".format(platform.release(),
                                                 platform.version(), sys_prop['Processor']))
    else:
        print("Unknown OS ({0})".format(sys_prop['Processor']))

    print("\n")
# required libraries
    print("###########")
    print(" Required Libraries ")
    print("###########")

    for sys_info in ['Python', 'NumPy', 'SciPy',
              'matplotlib', 'Astropy', 'Pandas']:
        print('{0}: {1}'.format(sys_info, sys_prop[sys_info]))

    print("\n")

# recommended
    print("###########")
    print(" Recommended Libraries ")
    print("###########")

    for sys_info in ['beautifulsoup', 'PyQt', 'SUDS',
                     'Sqlalchemy', 'Requests']:
        print('{0}: {1}'.format(sys_info, sys_prop[sys_info]))
