import pytest

import astropy.units as u
from astropy.table import Column, MaskedColumn, Table
from astropy.tests.helper import assert_quantity_allclose

from sunpy.data.test import get_test_filepath
from sunpy.io.special import srs

filenames = [
    # Mixedcase files
    {'file': '20150906SRS.txt', 'rows': 5},
    {'file': '20150306SRS.txt', 'rows': 4},
    {'file': '20150101SRS.txt', 'rows': 9},
    {'file': '20100621SRS.txt', 'rows': 3},  # This is a corrected copy
    # Uppercase files
    {'file': '19960106SRS.txt', 'rows': 4},  # inc. spurious `NNN` on final line
    {'file': '19960430SRS.txt', 'rows': 1},
    {'file': '19960513SRS.txt', 'rows': 4},  # inc. empty `COMMENT` column
    {'file': '20000922SRS.txt', 'rows': 10},  # inc. line with `III.`
    {'file': '20000927SRS.txt', 'rows': 9},  # inc. line with `EFFECTIVE 2 OCT`
    {'file': '20001001SRS.txt', 'rows': 12},  # inc. line with `COMMENT`
    {'file': '20020624SRS.txt', 'rows': 14},  # inc. line with `PLAIN`
    {'file': '20020628SRS.txt', 'rows': 13},  # inc. line with `This message`
]

COORDINATES = [{'text': 'N10W05', 'latitude': 10, 'longitude': 5},
               {'text': 'N89E00', 'latitude': 89, 'longitude': 0},
               {'text': 'S33E02', 'latitude': -33, 'longitude': -2},
               {'text': 'S01', 'latitude': -1, 'longitude': None}]

LOCATION = Column(data=[x['text'] for x in COORDINATES], name='Location')
LONGLAT = Table()
LONGLAT.add_column(MaskedColumn(data=[x['longitude'] for x in COORDINATES], name='Longitude',
                                unit=u.deg, mask=True, fill_value=-99999))
LONGLAT.add_column(MaskedColumn(data=[x['latitude'] for x in COORDINATES], name='Latitude',
                                unit=u.deg, fill_value=-99999))


@pytest.mark.filterwarnings('ignore:dropping mask in Quantity column')
@pytest.mark.parametrize(
    ('path', 'number_of_rows'),
    [(get_test_filepath(elem['file']), elem['rows']) for elem in filenames]
)
def test_number_of_rows(path, number_of_rows):
    table = srs.read_srs(path)
    assert len(table) == number_of_rows


@pytest.mark.parametrize(('text', 'longitude'),
                         [(elem['text'], elem['longitude']) for elem in COORDINATES])
def test_parse_longitude(text, longitude):
    assert srs.parse_longitude(text) == longitude


@pytest.mark.parametrize(('text', 'latitude'),
                         [(elem['text'], elem['latitude']) for elem in COORDINATES])
def test_parse_latitude(text, latitude):
    assert srs.parse_latitude(text) == latitude


@pytest.mark.parametrize(('loc_column', 'exp_longitude', 'exp_latitude'),
                         [(LOCATION, LONGLAT['Longitude'], LONGLAT['Latitude'])])
def test_parse_location(loc_column, exp_longitude, exp_latitude):
    latitude, longitude = srs.parse_location(loc_column)
    assert_quantity_allclose(latitude, exp_latitude)
    assert_quantity_allclose(longitude, exp_longitude)
