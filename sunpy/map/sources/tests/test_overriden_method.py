import pytest
from sunpy.data.test import get_test_filepath
import sunpy.map
from sunpy.map.mapbase import SpatialPair
from sunpy.data.test import get_dummy_map_from_header, get_test_filepath

def test_eitmap_coordinate_system():
    eit_map = sunpy.map.Map(get_test_filepath("EIT/efz20040301.000010_s.fits"))
    assert eit_map.coordinate_system ==  SpatialPair(axis1='HPLN-TAN', axis2='HPLT-TAN')


def test_lasco_coordinate_system():
    lasco_map = get_dummy_map_from_header(get_test_filepath("lasco_c2_25299383_s.header"))
    assert lasco_map.coordinate_system ==  SpatialPair(axis1='HPLN-TAN', axis2='HPLT-TAN')

def test_xrt_coordinate_system():
    xrt_map = get_dummy_map_from_header(get_test_filepath("HinodeXRT.header"))
    assert xrt_map.coordinate_system ==  SpatialPair(axis1='HPLN-TAN', axis2='HPLT-TAN')

def test_trace_coordinate_system():
    trace_map = get_dummy_map_from_header(get_test_filepath("tsi20010130_025823_a2.header"))
    assert trace_map.coordinate_system ==  SpatialPair(axis1='HPLN-TAN', axis2='HPLT-TAN')

def test_sot_coordinate_system():
    sot_map = get_dummy_map_from_header(get_test_filepath("HinodeSOT.header"))
    assert sot_map.coordinate_system ==  SpatialPair(axis1='HPLN-TAN', axis2='HPLT-TAN')

