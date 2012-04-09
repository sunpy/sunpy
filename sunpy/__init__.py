"""
SunPy
=====

An open-source Python library for Solar Physics data analysis.

Classes
-------
Map
    A spatially-aware 2d data array
MapCube
    A spatially-aware 3d data array


Available subpackages
---------------------
map
    Methods relating to the solar maps

"""
from __future__ import absolute_import

__version__ = 0.1

import sunpy.map
from sunpy.map import make_map
from sunpy.map import Map
from sunpy.map.header import MapHeader
from sunpy.map.mapcube import MapCube
from sunpy.map.compositemap import CompositeMap
from sunpy.cm import *

import sys, os, tempfile
from configobj import ConfigObj

# Sample data
from sunpy.data.sample import (AIA_171_IMAGE, RHESSI_IMAGE, EIT_195_IMAGE, 
                               RHESSI_EVENT_LIST)

def _is_writable_dir(p):
    """
    p is a string pointing to a putative writable dir -- return True p
    is such a string, else False
    """
    try: p + ''  # test is string like
    except TypeError: return False
    try:
        t = tempfile.TemporaryFile(dir=p)
        t.write('1')
        t.close()
    except OSError: return False
    else: return True

def _get_home():
    """Find user's home directory if possible.
    Otherwise raise error.

    """
    path=''
    try:
        path=os.path.expanduser("~")
    except:
        pass
    if not os.path.isdir(path):
        for evar in ('HOME', 'USERPROFILE', 'TMP'):
            try:
                path = os.environ[evar]
                if os.path.isdir(path):
                    break
            except: pass
    if path:
        return path
    else:
        raise RuntimeError('please define environment variable $HOME')

#get_home = verbose.wrap('$HOME=%s', _get_home, always=False)

def get_configdir():
    """
    Return the string representing the configuration dir.
    The default is "HOME/.sunpy".  You can override this with the
    SUNPY_CONFIGDIR environment variable
    """

    configdir = os.environ.get('SUNPY_CONFIGDIR')
    if configdir is not None:
        if not _is_writable_dir(configdir):
            raise RuntimeError('Could not write to SUNPY_CONFIGDIR="%s"'%configdir)
        return configdir

    h = _get_home()
    p = os.path.join(_get_home(), '.sunpy')

    if os.path.exists(p):
        if not _is_writable_dir(p):
            raise RuntimeError("'%s' is not a writable dir; you must set %s/.sunpy to be a writable dir.  You can also set environment variable SUNPY_CONFIGDIR to any writable directory where you want matplotlib data stored "% (h, h))
    else:
        if not _is_writable_dir(h):
            raise RuntimeError("Failed to create %s/.sunpy; consider setting SUNPY_CONFIGDIR to a writable directory for sunpy configuration data"%h)

        os.mkdir(p)

    return p

#get_configdir = verbose.wrap('CONFIGDIR=%s', _get_configdir, always=False)

def read_configfile():
    """
    Read the sunpyrc configuration file. If one does not exists in the user's
    home directory then read in the defaults from module
    """

    config_filename = 'sunpyrc'
    config_path = get_configdir()
    
    # check if to see if defaults have been customized 
    # if not read in the defaults from the module
    if os.path.exists(config_path + '/' + config_filename):    
        config = ConfigObj(config_path + '/' + config_filename)
    else:
        module_dir = os.path.dirname(sunpy.__file__)
        config = ConfigObj(module_dir + '/data/sunpyrc')
        
    return config
   
configs = read_configfile()