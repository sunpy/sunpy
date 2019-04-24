# -*- coding: utf-8 -*-

import logging
import os.path

from astropy.logger import AstropyLogger

from sunpy import config, log

level_to_numeric = {'CRITICAL': 50, 'ERROR':40, 'WARNING': 30, 'INFO': 20, 'DEBUG': 10, 'NOTSET': 0}


def test_logger_name():
    assert log.name == 'sunpy'


def test_is_the_logger_there():
    assert isinstance(log, logging.Logger)
    assert isinstance(log, AstropyLogger)


def test_is_level_configed():
    """
    Test to make sure that the logger follows the config:

    log_level
    """
    config_level_numeric = level_to_numeric.get(config.get('logger', 'log_level'))
    assert log.getEffectiveLevel() == config_level_numeric


def test_is_log_to_file_configed():
    """
    Test to make sure that the logger follows the config:

    log_to_file, log_file_level, log_file_path
    """

    if config.get('logger', 'log_to_file') == 'True':
        #  there must be two handlers, one streaming and one to file.
        assert len(log.handlers) == 2
        #  one of the handlers must be FileHandler
        assert isinstance(log.handlers[0], logging.FileHandler) or isinstance(log.handlers[1], logging.FileHandler)
        fh = None
        if isinstance(log.handlers[0], logging.FileHandler):
            fh = log.handlers[0]

        if isinstance(log.handlers[1], logging.FileHandler):
            fh = log.handlers[1]

        if fh is not None:
            log_file_level = config.get('logger', 'log_file_level')
            assert level_to_numeric.get(log_file_level) == fh.level

            log_file_path = config.get('logger', 'log_file_path')
            assert os.path.basename(fh.baseFilename) == os.path.basename(log_file_path)


def test_origin():
    with log.log_to_list() as log_list:
        log.info('test1')

    assert log_list[0].origin == 'sunpy.util.tests.test_logger'
    assert log_list[0].message.startswith('test1')


def send_to_log(message, kind='INFO'):
    """
    A simple function to demonstrate the logger generating an origin.
    """
    if kind.lower() == 'info':
        log.info(message)
    elif kind.lower() == 'debug':
        log.debug(message)

# no obvious way to do the following
# TODO: test for the following configs  use_color, log_warnings, log_exceptions, log_file_format
