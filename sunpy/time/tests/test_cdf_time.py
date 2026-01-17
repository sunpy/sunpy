
import sys
from unittest import mock
import pytest
import numpy as np
from astropy.time import Time
from sunpy.time import parse_time

def test_parse_time_cdf_tt2000_scalar():
    pytest.importorskip("cdflib")
    # Example scalar TT2000
    # J2000 = 2000-01-01 12:00:00 TT
    tt2000 = 0
    t = parse_time(tt2000, format='cdf_tt2000')
    assert isinstance(t, Time)
    # Depending on version, format might be 'cdf_tt2000' or converted to 'isot' etc
    # cdflib 1.3.8 seems to return Time object with format='cdf_tt2000'
    # But checks might vary.
    assert t.shape == ()
    
def test_parse_time_cdf_tt2000_array():
    pytest.importorskip("cdflib")
    tt2000 = np.array([0, 1000000000]) # 0 and 1 sec
    t = parse_time(tt2000, format='cdf_tt2000')
    assert isinstance(t, Time)
    assert len(t) == 2
    # Verify values logic if possible, but mainly verification of integration.

def test_parse_time_cdf_epoch():
    pytest.importorskip("cdflib")
    # CDF_EPOCH: ms since 0000-01-01
    epoch = 63830592000000.0
    t = parse_time(epoch, format='cdf_epoch')
    assert isinstance(t, Time)

def test_missing_cdflib_error():
    # We simulate missing cdflib by patching sys.modules
    with mock.patch.dict(sys.modules):
        # Determine all cdflib modules
        modules_to_remove = [m for m in sys.modules if m.startswith("cdflib")]
        for m in modules_to_remove:
            del sys.modules[m]
        
        # Also need to make sure import fails. 
        # assigning None to sys.modules['cdflib'] usually causes ImportError or ModuleNotFoundError
        sys.modules['cdflib'] = None
        sys.modules['cdflib.epochs_astropy'] = None
        
        with pytest.raises(ImportError, match="cdflib must be installed"):
             parse_time(0, format='cdf_tt2000')

