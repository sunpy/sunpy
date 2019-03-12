import logging

from astropy.logger import AstropyLogger


def _init_log(config=None):
    """Initializes the SunPy log--in most circumstances this is called
    automatically when importing sunpy.
    """

    global log

    orig_logger_cls = logging.getLoggerClass()
    logging.setLoggerClass(AstropyLogger)
    try:
        log = logging.getLogger('sunpy')
        if config is not None:
            conf = _config_to_loggerConf(config)
        log._set_defaults()
    finally:
        logging.setLoggerClass(orig_logger_cls)

    return log


def _config_to_loggerConf(config):
    """Set the log configuration based on the user-provided config."""

    if config.has_section('logger'):
        from astropy.logger import Conf as LoggerConf
        conf = LoggerConf()
        option_list = ['log_level', 'use_color', 'log_warnings', 'log_exceptions', 'log_to_file', 'log_file_path',
                       'log_file_level', 'log_file_format']
        for this_option in option_list:
            if config.has_option('logger', this_option):
                setattr(conf, this_option, config.get('logger', this_option))
    return conf
