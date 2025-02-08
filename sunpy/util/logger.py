import os
import sys
import logging

from astropy.logger import AstropyLogger

from sunpy.util.exceptions import SunpyWarning


class SunpyLogger(AstropyLogger):
    """
    This class is used to set up logging in SunPy.

    This inherits the logging enhancements of `~astropy.logger.AstropyLogger`.
    This logger is able to capture and log warnings that are based on
    `~sunpy.util.exceptions.SunpyWarning`. Other warnings will be ignored and
    passed on to other loggers (e.g., from Astropy).
    """
    # Override the existing _showwarning() to capture SunpyWarning instead of AstropyWarning

    def _showwarning(self, *args, **kwargs):

        # Bail out if we are not catching a warning from SunPy
        if not isinstance(args[0], SunpyWarning):
            return self._showwarning_orig(*args, **kwargs)

        warning = args[0]
        # Deliberately not using isinstance here: We want to display
        # the class name only when it's not the default class,
        # SunpyWarning. The name of subclasses of SunpyWarning should
        # be displayed.
        if type(warning) not in (SunpyWarning,):
            message = f'{warning.__class__.__name__}: {args[0]}'
        else:
            message = str(args[0])

        mod_path = args[2]
        # Now that we have the module's path, we look through sys.modules to
        # find the module object and thus the fully-package-specified module
        # name. The module.__file__ is the original source file name.
        mod_name = None
        mod_path, ext = os.path.splitext(mod_path)
        for name, mod in list(sys.modules.items()):
            try:
                # Believe it or not this can fail in some cases:
                # https://github.com/astropy/astropy/issues/2671
                path = os.path.splitext(getattr(mod, '__file__', ''))[0]
            except Exception:
                continue
            if path == mod_path:
                mod_name = mod.__name__
                break

        if mod_name is not None:
            self.warning(message, extra={'origin': mod_name})
        else:
            self.warning(message)


def _init_log(config=None):
    """
    Initializes the SunPy log.

    In most circumstances this is called automatically when importing
    SunPy. This code is based on that provided by Astropy see
    "licenses/ASTROPY.rst".
    """
    orig_logger_cls = logging.getLoggerClass()
    logging.setLoggerClass(SunpyLogger)
    try:
        log = logging.getLogger('sunpy')
        if config is not None:
            _config_to_loggerConf(config)
        log._set_defaults()
    finally:
        logging.setLoggerClass(orig_logger_cls)

    return log


def _config_to_loggerConf(config):
    """
    Translates a user-provided config to ~`astropy.logger.LoggerConf`.
    """

    if config.has_section('logger'):
        from astropy.logger import Conf as LoggerConf
        conf = LoggerConf()
        loggerconf_option_list = ['log_level', 'use_color', 'log_warnings', 'log_exceptions', 'log_to_file',
                                  'log_file_path', 'log_file_level', 'log_file_format']
        for this_option in loggerconf_option_list:
            if config.has_option('logger', this_option):
                setattr(conf, this_option, config.get('logger', this_option))
    return conf
