"""
This module implements tests for SRS Reader.
"""
import os

import pytest
import numpy as np
import astropy.units as u
from astropy.table import Table, Column, MaskedColumn
from astropy.tests.helper import assert_quantity_allclose

from sunpy.io.special import srs
import sunpy.data.test

testpath = sunpy.data.test.rootdir

filenames = [{'file': '20150906SRS.txt', 'rows': 5},
             {'file': '20150306SRS.txt', 'rows': 4},
             {'file': '20150101SRS.txt', 'rows': 9}]

COORDINATES = [{'text': 'N10W05', 'latitude': 10,  'longitude': 5},
               {'text': 'N89E00', 'latitude': 89,  'longitude': 0},
               {'text': 'S33E02', 'latitude': -33, 'longitude': -2},
               {'text': 'S01', 'latitude': -1, 'longitude': None}]

LOCATION = Column(data=[x['text'] for x in COORDINATES], name='Location')
LONGLAT = Table()
LONGLAT.add_column(MaskedColumn(data=[x['longitude'] for x in COORDINATES], name='Longitude',
                                unit=u.deg, mask=True))
LONGLAT.add_column(MaskedColumn(data=[x['latitude'] for x in COORDINATES], name='Latitude',
                                unit=u.deg))


@pytest.mark.parametrize("path, number_of_rows",
                         [(os.path.join(testpath, elem['file']), elem['rows'])
                          for elem in filenames])
def test_number_of_rows(path, number_of_rows):
    table = srs.read_srs(path)
    assert len(table) == number_of_rows


@pytest.mark.parametrize("text, longitude",
                         [(elem['text'], elem['longitude']) for elem in COORDINATES])
def test_parse_longitude(text, longitude):
    assert srs.parse_longitude(text) == longitude


@pytest.mark.parametrize("text, latitude",
                         [(elem['text'], elem['latitude']) for elem in COORDINATES])
def test_parse_latitude(text, latitude):
    assert srs.parse_latitude(text) == latitude


@pytest.mark.parametrize("loc_column, exp_longitude, exp_latitude",
                         [(LOCATION, LONGLAT['Longitude'], LONGLAT['Latitude'])])
def test_parse_location(loc_column, exp_longitude, exp_latitude):
    latitude, longitude = srs.parse_location(loc_column)
    assert_quantity_allclose(latitude, exp_latitude)
    assert_quantity_allclose(longitude, exp_longitude)
