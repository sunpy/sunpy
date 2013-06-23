from __future__ import absolute_import

import platform
import datetime

import sunpy

__all__ = ['system_info']

def system_info():
    """Prints system information.
    
    Prints information about the runtime environment that SunPy lives in.
    Information about the OS, architecture, Python, and all major dependencies
    are included.
    
    The goal of this function is to provide enough information for someone
    running SunPy code or replicating a bug to setup a comparible environment
    to that which was originally used.
    
    Author: `Keith Hughitt <keith.hughitt@nasa.gov>` 
    """
    
    print("==========================================================")
    print(" SunPy Installation Information\n")
    print(" " + datetime.datetime.utcnow().strftime("%A, %d. %B %Y %I:%M%p UT"))
    print("==========================================================\n")
    
    system = platform.system()
    proc = platform.processor()
    
    print("###########")
    print(" General")
    print("###########")
    
    # OS and architecture information
    if system == "Linux":
        distro = " ".join(platform.linux_distribution())
        print("OS: %s (Linux %s %s)" %  (distro, platform.release(), proc))
    elif system == "Darwin":
        print("OS: Mac OS X %s (%s)" %  (platform.mac_ver()[0], proc))
    elif system == "Windows":
        print("OS: Windows %s %s (%s)" %  (platform.release(), 
                                        platform.version(), proc))
    else:
        print ("Unknown OS (%s)" % proc)
        
    # Python version
    arch = platform.architecture()[0]
    print("Python: %s (%s)\n" % (platform.python_version(), arch))
    
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
        from pyfits import __version__ as pyfits_version
    except ImportError:
        pyfits_version = "NOT INSTALLED"
        
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

    print("####################")
    print(" Required libraries")
    print("####################")
    
    print("SunPy: %s" % sunpy.__version__)

    print("NumPy: %s" % numpy_version)
    print("SciPy: %s" % scipy_version)
    print("Matplotlib: %s" % matplotlib_version)
    print("PyFITS: %s" % pyfits_version)
    print("pandas: %s" % pandas_version)

    print("")

    print("#######################")
    print(" Recommended libraries")
    print("#######################")

    print("beautifulsoup4: %s" % bs4_version)
    print("PyQt: %s" % pyqt_version)
    print("SUDS: %s" % suds_version)

    print("")
