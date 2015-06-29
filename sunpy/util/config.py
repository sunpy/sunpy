"""SunPy configuration file functionality"""
import os
import tempfile
import ConfigParser

import sunpy as spy

__all__ = ['load_config', 'print_config']

def load_config():
    """
    Read the sunpyrc configuration file. If one does not exists in the user's
    home directory then read in the defaults from module
    """
    config = ConfigParser.SafeConfigParser()

    # Get locations of SunPy configuration files to be loaded
    config_files = _find_config_files()

    # Read in configuration files
    config.read(config_files)

    # Specify the working directory as a default so that the user's home
    # directory can be located in an OS-independent manner
    if not config.has_option('general', 'working_dir'):
        config.set('general', 'working_dir', os.path.join(_get_home(), "sunpy"))

    # Specify the database url as a default so that the user's home
    # directory can be located in an OS-independent manner
    if not config.has_option('database', 'url'):
        config.set('database', 'url', "sqlite:///" + os.path.join(_get_home(), "sunpy/sunpydb.sqlite"))

    # Use absolute filepaths and adjust OS-dependent paths as needed
    filepaths = [
        ('downloads', 'download_dir'),
        ('downloads', 'sample_dir')
    ]
    _fix_filepaths(config, filepaths)

    # check for sunpy working directory and create it if it doesn't exist
    if not os.path.isdir(config.get('downloads', 'download_dir')):
        os.mkdir(config.get('downloads', 'download_dir'))

    if not os.path.isdir(config.get('downloads', 'sample_dir')):
        os.mkdir(config.get('downloads', 'sample_dir'))

    return config

def print_config():
    """Print current configuration options"""
    print("FILES USED:")
    for file_ in _find_config_files():
        print("  " + file_)

    print ("\nCONFIGURATION:")
    for section in spy.config.sections():
        print("  [{0}]".format(section))
        for option in spy.config.options(section):
            print("  {0} = {1}".format(option, spy.config.get(section, option)))
        print("")

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

def _find_config_files():
    """Finds locations of SunPy configuration files"""
    config_files = []
    config_filename = 'sunpyrc'

    # find default configuration file
    module_dir = os.path.dirname(spy.__file__)
    config_files.append(os.path.join(module_dir, 'data', 'sunpyrc'))

    # if a user configuration file exists, add that to list of files to read
    # so that any values set there will override ones specified in the default
    # config file
    config_path = _get_user_configdir()

    if os.path.exists(os.path.join(config_path, config_filename)):
        config_files.append(os.path.join(config_path, config_filename))

    return config_files

def _get_user_configdir():
    """
    Return the string representing the configuration dir.
    The default is "HOME/.sunpy".  You can override this with the
    SUNPY_CONFIGDIR environment variable
    """
    configdir = os.environ.get('SUNPY_CONFIGDIR')

    if configdir is not None:
        if not _is_writable_dir(configdir):
            raise RuntimeError('Could not write to SUNPY_CONFIGDIR="{0}"'.format(
                               configdir))
        return configdir

    h = _get_home()
    p = os.path.join(_get_home(), '.sunpy')

    if os.path.exists(p):
        if not _is_writable_dir(p):
            raise RuntimeError("'{0}' is not a writable dir; you must set {1}/."
                               "sunpy to be a writable dir.  You can also set "
                               "environment variable SUNPY_CONFIGDIR to any "
                               "writable directory where you want matplotlib "
                               "data stored ".format(h, h))
    else:
        if not _is_writable_dir(h):
            raise RuntimeError("Failed to create {0}/.sunpy; consider setting "
                               "SUNPY_CONFIGDIR to a writable directory for "
                               "sunpy configuration data".format(h))

        os.mkdir(p)

    return p

def _fix_filepaths(config, filepaths):
    """Converts relative filepaths to absolute filepaths"""
    # Parse working_dir
    working_dir = _expand_filepath(config.get("general", "working_dir"))
    config.set('general', 'working_dir', working_dir)

    for f in filepaths:
        val = config.get(*f)

        filepath = _expand_filepath(val, working_dir)

        # Create dir if it doesn't already exist
        if not os.path.isdir(filepath):
            os.makedirs(filepath)

        # Replace config value with full filepath
        params = f + (filepath,)
        config.set(*params)

def _expand_filepath(filepath, working_dir=""):
    """Checks a filepath and expands it if necessary"""
    # Expand home directory
    if filepath[0] == "~":
        return os.path.abspath(os.path.expanduser(filepath))
    # Check for /tmp
    elif filepath == "/tmp":
        return tempfile.gettempdir()
    # Relative filepaths
    elif not filepath.startswith("/"):
        return os.path.join(working_dir, filepath)
    # Absolute filepath
    else:
        return filepath
