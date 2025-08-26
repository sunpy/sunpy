
import logging
import os.path
import warnings

import pytest

from astropy.utils.exceptions import AstropyUserWarning

from sunpy import config, log
from sunpy.util.exceptions import SunpyUserWarning
from sunpy.util.logger import SunpyLogger

level_to_numeric = {'CRITICAL': 50, 'ERROR': 40,
                    'WARNING': 30, 'INFO': 20, 'DEBUG': 10, 'NOTSET': 0}


def test_logger_name():
    assert log.name == 'sunpy'


def test_is_the_logger_there():
    assert isinstance(log, logging.Logger)
    assert isinstance(log, SunpyLogger)


def test_is_level_configured():
    """
    Test to make sure that the logger follows the config:

    log_level
    """
    config_level_numeric = level_to_numeric.get(config.get('logger', 'log_level'))
    assert log.getEffectiveLevel() == config_level_numeric


def test_is_log_to_file_configured():
    """
    Test to make sure that the logger follows the config:

    log_to_file, log_file_level, log_file_path
    """

    if config.get('logger', 'log_to_file') == 'True':
        #  there must be two handlers, one streaming and one to file.
        assert len(log.handlers) == 2
        #  one of the handlers must be FileHandler
        assert isinstance(log.handlers[0], logging.FileHandler) or isinstance(
            log.handlers[1], logging.FileHandler)
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


# Most of the logging functionality is tested in Astropy's tests for AstropyLogger

def test_sunpy_warnings_logging():
    # Test that our logger intercepts our warnings but not Astropy warnings

    # First, disable our warnings logging
    # We need to do this manually because pytest has overwritten warnings.showwarning()
    log._showwarning_orig, previous = None, log._showwarning_orig

    # Without warnings logging
    with pytest.warns(SunpyUserWarning, match="This warning should not be captured") as warn_list:
        with log.log_to_list() as log_list:
            warnings.warn("This warning should not be captured", SunpyUserWarning)
    assert len(log_list) == 0
    assert len(warn_list) == 1

    # With warnings logging, making sure that Astropy warnings are not intercepted
    with pytest.warns(AstropyUserWarning, match="This warning should not be captured") as warn_list:  # NOQA: PT031
        log.enable_warnings_logging()
        with log.log_to_list() as log_list:
            warnings.warn("This warning should be captured", SunpyUserWarning)
            warnings.warn("This warning should not be captured", AstropyUserWarning)
        log.disable_warnings_logging()
    assert len(log_list) == 1
    assert len(warn_list) == 1
    assert log_list[0].levelname == 'WARNING'
    assert log_list[0].message.startswith("SunpyUserWarning: This warning should be captured")
    assert log_list[0].origin == "sunpy.util.tests.test_logger"

    # Restore the state of warnings logging prior to this test
    log._showwarning_orig = previous
