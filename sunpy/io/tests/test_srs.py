import pytest

import astropy.units as u
from astropy.table import Column, MaskedColumn, Table
from astropy.tests.helper import assert_quantity_allclose

from sunpy.data.test import get_test_filepath
from sunpy.io.special import srs

filenames = [{'file': 'SRS/20150906SRS.txt', 'rows': 5},
             {'file': 'SRS/20150306SRS.txt', 'rows': 4},
             {'file': 'SRS/20150101SRS.txt', 'rows': 9},
             {'file': 'SRS/20100621SRS.txt', 'rows': 3},  # This is a corrected copy
             # Uppercase files
             {'file': 'SRS/19960106SRS.txt', 'rows': 4},  # inc. spurious `NNN` on final line
             {'file': 'SRS/19960430SRS.txt', 'rows': 1},
             {'file': 'SRS/19960513SRS.txt', 'rows': 4},  # inc. empty `COMMENT` column
             {'file': 'SRS/20000922SRS.txt', 'rows': 10},  # inc. line with `III.`
             {'file': 'SRS/20000927SRS.txt', 'rows': 9},  # inc. line with `EFFECTIVE 2 OCT`
             {'file': 'SRS/20001001SRS.txt', 'rows': 12},  # inc. line with `COMMENT`
             {'file': 'SRS/20020624SRS.txt', 'rows': 14},  # inc. line with `PLAIN`
             {'file': 'SRS/20020628SRS.txt', 'rows': 13}]  # inc. line with `This message`

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

# test _try_drop_empty_column (tdec)


@pytest.fixture
def data_lines():
    return ['NMBR  LOCATION  LO  COMMENT', '8000  S14W96   232']


@pytest.fixture
def expected_pattern_dict():
    return {
        'Nmbr': r'^\d+$',
        'Location': r'^(?:[NESW](?:\d{2})){1,2}$',
        'Lo': r'^\d+$',
        'Comment': r'^[a-zA-Z]+$',
        'Additional Number': r'^\d+$',
    }


@pytest.fixture
def column_name_to_drop():
    return 'COMMENT'


def test_tdec(data_lines, expected_pattern_dict, column_name_to_drop):
    result = srs._try_drop_empty_column(column_name_to_drop, data_lines, expected_pattern_dict)
    out = data_lines.copy()
    out[0] = ' '.join([col.title() for col in out[0].split()][:-1])
    assert result == out


def test_tdec_smallest_dict_example(data_lines, expected_pattern_dict, column_name_to_drop):
    """
    Smallest possible dictionary given the ``data_lines``
    """
    keys_to_keep = ['Nmbr', 'Location', 'Lo']
    filtered_dict = {key: value for key, value in expected_pattern_dict.items() if key in keys_to_keep}
    result = srs._try_drop_empty_column(column_name_to_drop, data_lines, filtered_dict)
    out = data_lines.copy()
    out[0] = ' '.join([col.title() for col in out[0].split()][:-1])
    assert result == out


def test_tdec_too_small_dict_example(data_lines, expected_pattern_dict, column_name_to_drop):
    """
    Columns aren't a subset of the dictionary: no pattern to match for column `LO`
    """
    keys_to_keep = ['Nmbr', 'Location']
    filtered_dict = {key: value for key, value in expected_pattern_dict.items() if key in keys_to_keep}
    with pytest.raises(ValueError, match="The remaining columns are not a subset of the columns in ``pattern_dict``."):
        _ = srs._try_drop_empty_column(column_name_to_drop, data_lines, filtered_dict)


def test_tdec_no_data(expected_pattern_dict, column_name_to_drop):
    """
    No data associated with header
    """
    data_lines = ['NMBR  LOCATION  LO  COMMENT', 'NONE']
    out = data_lines.copy()
    out[0] = ' '.join([col.title() for col in out[0].split()][:-1])
    assert out == srs._try_drop_empty_column(column_name_to_drop, data_lines, expected_pattern_dict)


def test_tdec_col_not_empty(expected_pattern_dict, column_name_to_drop):
    """
    Column not empty
    """
    data_lines = [
        'NMBR  LOCATION  LO  COMMENT',
        '8000  S14W96   232',
        '3020 N20E20  210 additional_info',
    ]
    with pytest.raises(ValueError, match="not all rows have the same number of values as the remaining columns."):
        _ = srs._try_drop_empty_column(column_name_to_drop, data_lines, expected_pattern_dict)


def test_tdec_col_not_match_pattern(expected_pattern_dict, column_name_to_drop):
    """
    `LOCATION` column doesn't match the correct pattern
    """
    data_lines = [
        'NMBR  LOCATION  LO  COMMENT',
        '8000  S14W926   232',
    ]
    with pytest.raises(ValueError, match="not all rows match the provided pattern."):
        _ = srs._try_drop_empty_column(column_name_to_drop, data_lines, expected_pattern_dict)


def test_tdec_colname_not_exist(data_lines, expected_pattern_dict):
    """
    Try remove column name that does not exist in data
    """
    with pytest.raises(ValueError, match="The column 'nonexistent_column_name' does not exist."):
        _ = srs._try_drop_empty_column('nonexistent_column_name', data_lines, expected_pattern_dict)
