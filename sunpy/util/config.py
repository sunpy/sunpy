"""SunPy configuration file functionality"""
import os
import sunpy
from configobj import ConfigObj

def _is_writable_dir(p):
    """Checks to see if a directory is writable"""
    return os.path.isdir(p) and os.access(p, os.W_OK)

def _get_home():
    """Find user's home directory if possible.
    Otherwise raise error.

    """
    path = path=os.path.expanduser("~")

    if not os.path.isdir(path):
        for evar in ('HOME', 'USERPROFILE', 'TMP'):
            try:
                path = os.environ[evar]
                if os.path.isdir(path):
                    break
            except KeyError:
                pass
    if path:
        return path
    else:
        raise RuntimeError('please define environment variable $HOME')

def _get_configdir():
    """
    Return the string representing the configuration dir.
    The default is "HOME/.sunpy".  You can override this with the
    SUNPY_CONFIGDIR environment variable
    """
    configdir = os.environ.get('SUNPY_CONFIGDIR')
    
    if configdir is not None:
        if not _is_writable_dir(configdir):
            raise RuntimeError('Could not write to SUNPY_CONFIGDIR="%s"' % 
                               configdir)
        return configdir

    h = _get_home()
    p = os.path.join(_get_home(), '.sunpy')

    if os.path.exists(p):
        if not _is_writable_dir(p):
            raise RuntimeError("'%s' is not a writable dir; you must set %s/."
                               "sunpy to be a writable dir.  You can also set "
                               "environment variable SUNPY_CONFIGDIR to any "
                               "writable directory where you want matplotlib "
                               "data stored " % (h, h))
    else:
        if not _is_writable_dir(h):
            raise RuntimeError("Failed to create %s/.sunpy; consider setting "
                               "SUNPY_CONFIGDIR to a writable directory for "
                               "sunpy configuration data" % h)

        os.mkdir(p)

    return p

def read_configfile():
    """
    Read the sunpyrc configuration file. If one does not exists in the user's
    home directory then read in the defaults from module
    """
    config_filename = 'sunpyrc'
    config_path = _get_configdir()
    
    # check if to see if defaults have been customized 
    # if not read in the defaults from the module
    if os.path.exists(config_path + '/' + config_filename):    
        config = ConfigObj(config_path + '/' + config_filename)
    else:
        module_dir = os.path.dirname(sunpy.__file__)
        config = ConfigObj(module_dir + '/data/sunpyrc')
        
    return config