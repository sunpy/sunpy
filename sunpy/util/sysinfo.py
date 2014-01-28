from __future__ import absolute_import

import platform
import datetime

import sunpy

__all__ = ['system_info', 'sys_prop_print']


    
def system_info():
    """
    

    """    
    
    
    system = platform.system()
    proc = platform.processor()
    
    
    try:
        from sunpy.version import version as sunpy_version
        from sunpy.version import git_description as sunpy_git_description
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
    except:
        suds_version = "NOT INSTALLED"
        
    sys_prop = {'Time' : datetime.datetime.utcnow().strftime("%A, %d. %B %Y %I:%M%p UT"),
                'System' : platform.system(), 'Processor' : platform.processor(), 
                'Arch' : platform.architecture()[0], "Python" : platform.python_version(),
                'SunPy': sunpy_version, 
                'SunPy_git' : sunpy_git_description, 'NumPy' : numpy_version, 
                'SciPy' : scipy_version, 'matplotlib' : matplotlib_version,
                'Astropy' : astropy_version, 'Pandas' : pandas_version, 
                'beautifulsoup' : bs4_version, 'PyQt' : pyqt_version,
                'SUDS' : suds_version
                }
    return sys_prop
    
def sys_prop_print():
    sys_prop = system_info()
    
# title
    print("==========================================================")
    print(" SunPy Installation Information\n")
    print("==========================================================\n")  
        
    
# general properties      
    print("###########")
    print(" General")
    print("###########")   
    # OS and architecture information
    
    for i in ['Time', 'System', 'Processor', 'Arch']:
        print i + ' : ' + sys_prop[i]
        
    if sys_prop['System'] == "Linux":
        distro = " ".join(platform.linux_distribution())
        print("OS: %s (Linux %s %s)" %  (distro, platform.release(), sys_prop['Processor']))
    elif sys_prop['System'] == "Darwin":
        print("OS: Mac OS X %s (%s)" %  (platform.mac_ver()[0], sys_prop['Processor']))
    elif sys_prop['System'] == "Windows":
        print("OS: Windows %s %s (%s)" %  (platform.release(), 
                                        platform.version(), sys_prop['Processor']))
    else:
        print ("Unknown OS (%s)" % sys_prop['Processor'])

    
    print "\n"
# required libraries
    print("###########")
    print(" Required Libraries ")
    print("###########")
    
    for i in ['Python', 'SunPy', 'SunPy_git', 'NumPy', 'SciPy',
              'matplotlib', 'Astropy', 'Pandas']:
        print i + ' : ' + sys_prop[i]
    
    print "\n"
    
# recommended
    print("###########")    
    print(" Recommended Libraries ") 
    print("###########")
        
    for i in ['beautifulsoup', 'PyQt', 'SUDS']:
        print i + ' : ' + sys_prop[i]