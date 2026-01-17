import pytest

from astropy.time import Time

from sunpy.time import parse_time

# Critical: Skip entire file if cdflib is missing
cdflib = pytest.importorskip("cdflib")
# Also skip if cdflib is present but too old to have epochs_astropy
try:
    import cdflib.epochs_astropy
except ImportError:
    pytest.skip("cdflib installed but too old (missing epochs_astropy)", allow_module_level=True)

def test_parse_time_cdf_tt2000_scalar():
    t_ns = 599572869184000000
    t = parse_time(t_ns, format='cdf_tt2000')
    assert isinstance(t, Time)
    assert t.format == 'cdf_tt2000'
    assert t.value == t_ns
    assert t.scale == 'tt'

def test_parse_time_cdf_tt2000_array():
    t_ns_array = [599572869184000000, 599572870184000000]
    t = parse_time(t_ns_array, format='cdf_tt2000')
    assert isinstance(t, Time)
    assert len(t) == 2
    assert t[0].value == t_ns_array[0]
    assert t[1].value == t_ns_array[1]

def test_parse_time_cdf_epoch():
    t_epoch = 63806836800000.0
    t = parse_time(t_epoch, format='cdf_epoch')
    assert isinstance(t, Time)
    assert t.format == 'cdf_epoch'
    assert t.value == t_epoch

import sys
from unittest import mock


def test_missing_cdflib_error():
    # We must patch sys.modules to simulate cdflib being missing.
    # The key is to ensure we also remove it from the conversion function's scope if it was already imported.
    with mock.patch.dict(sys.modules, {"cdflib": None, "cdflib.epochs_astropy": None}):
        with pytest.raises(ImportError, match="cdflib must be installed"):
            parse_time(0, format='cdf_tt2000')
