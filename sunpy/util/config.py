"""SunPy configuration file functionality"""
import os
import sunpy
import tempfile
import ConfigParser

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

def _fix_filepaths(config):
    """Converts relative filepaths to absolute filepaths"""
    # Filepath config parameters
    filepaths = [('downloads', 'download_dir')]

    for f in filepaths:
        val = config.get(*f)
        
        # Check for /tmp
        if val == "/tmp":
            val = tempfile.gettempdir()

        # Expand filepaths
        params = f + (os.path.abspath(os.path.expanduser(val)),)
        config.set(*params)

def read_configfile():
    """
    Read the sunpyrc configuration file. If one does not exists in the user's
    home directory then read in the defaults from module
    """
    config = ConfigParser.ConfigParser()
    
    # determine location of config file
    config_filename = 'sunpyrc'
    config_path = _get_configdir()
    
    # check if to see if defaults have been customized 
    # if not read in the defaults from the module
    if os.path.exists(config_path + '/' + config_filename):
        filepath = config_path + '/' + config_filename
    else:
        module_dir = os.path.dirname(sunpy.__file__)
        filepath = module_dir + '/data/sunpyrc'
    
    # Read in configuration
    config.readfp(open(filepath))
    
    # Use absolute filepaths and adjust OS-dependent paths as needed
    _fix_filepaths(config)
        
    return config

if __name__ == "__main__":
    import sunpy
    sunpy.util.system_info()