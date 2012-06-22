#!/usr/bin/env python
"""
Mac OS X SunPy installation script

This helper script checks the version of Python currently available, and if
sufficient installs SunPy and it's prereqs.
"""
import os
import sys

# Check python version
if sys.version_info < (2, 6):
    print("Please install the latest 2.x version of Python from www.python.org "
          "before continuing.")
    sys.exit()

# Check for GCC (Xcode)
if os.system('which gcc') is not 0:
    print("Please install XCode before continuing.")
    sys.exit()
    
# Install and update Homebrew
os.system('/usr/bin/ruby -e "$(/usr/bin/curl -fsSL https://raw.github.com/mxcl/homebrew/master/Library/Contributions/install_homebrew.rb)"')
os.system('brew doctor')

# Install pip
os.system('sudo easy_install pip')

# Install main compilation reqs
os.system('brew -v  install gfortran pkgconfig git openjpeg qt')

# Install rest of SunPy dependencies
os.system('sudo pip install --upgrade distribute')
os.system('sudo pip install --upgrade ipython')
os.system('sudo pip install --upgrade numpy')
os.system('sudo pip install --upgrade scipy')
os.system('sudo pip install --upgrade pyfits')
os.system('sudo pip install --upgrade suds')
os.system('sudo pip install --upgrade pandas')
os.system('sudo pip install --upgrade matplotlib')

# Install SunPy
os.system('sudo pip install git+git://github.com/sunpy/sunpy.git')
