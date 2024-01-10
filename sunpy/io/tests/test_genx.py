import datetime

import numpy as np
import pytest

from sunpy.data.test import get_test_filepath
from sunpy.io.special import genx

TESTING = genx.read_genx(get_test_filepath('generated_sample.genx'))


def test_skeleton():
    # Top level
    toplevel_dims = {'MYTEXT': 63, 'MYTEXT_ARRAY': 3, 'MYTEXT_ARRAY_DIMENSION': (2, 3),
                     'MYNUMBER': 1, 'MYNUMBER_ARRAY': 3, 'MYNUMBER_ARRAY_DIMENSION': (2, 3, 4, 5),
                     'MYUINT': 1, 'MYSTRUCTURE': 14,  # the elements inside the OrderedDict
                     'MYSTRUCTURE_ARRAY': 6, 'HEADER': 5}
    assert sorted(list(TESTING.keys())) == sorted(list(toplevel_dims.keys()))
    for key, val in toplevel_dims.items():
        if isinstance(val, tuple):
            assert TESTING[key].shape == tuple(reversed(val))
        else:
            if val > 1:
                assert len(TESTING[key]) == val
            else:
                assert isinstance(TESTING[key], int)


def test_array_elements_values():
    np.testing.assert_allclose(TESTING['MYSTRUCTURE']['MYFARRAY'], np.arange(3.))
    np.testing.assert_allclose(TESTING['MYSTRUCTURE']['MYFARRAYD'][:, 0], np.arange(6., step=2))
    assert TESTING['MYSTRUCTURE']['MYDARRAYD'][1, 2] == 5.
    assert TESTING['MYSTRUCTURE']['NESTEDSTRUCT']['MYLARRAYD'][3, 0, 1] == 19
    np.testing.assert_allclose(TESTING['MYSTRUCTURE']['NESTEDSTRUCT']
                               ['MYLARRAYD'][2, :, 0], np.arange(12, 17, step=2))
    assert TESTING['MYSTRUCTURE']['MYCARRAY'][1] == complex(1, -9)
    assert TESTING['MYSTRUCTURE']['MYDCARRAY'][2] == complex(12, 1)
    assert TESTING['MYSTRUCTURE']['NESTEDSTRUCT']['MYUL64NUMBER'] == 18446744073709551615
    assert TESTING['MYSTRUCTURE']['NESTEDSTRUCT']['MYL64NUMBER'] == 9223372036854775807


@pytest.mark.parametrize(('slice', 'value'), [((0, 0, 0, 0), 0),
                                              ((4, 0, 0, 0), 96),
                                              ((0, 2, 2, 0), 16),
                                              ((0, 3, 2, 0), 22),
                                              ((4, 3, 2, 0), 118)])
def test_value_slice(slice, value):
    assert TESTING['MYNUMBER_ARRAY_DIMENSION'][slice] == value


@pytest.mark.parametrize(('myarray', 'dtype'), [(TESTING['MYNUMBER_ARRAY'], np.int16),
                                                (TESTING['MYNUMBER_ARRAY_DIMENSION'], np.int16),
                                                (TESTING['MYSTRUCTURE']['MYFARRAY'], np.float32),
                                                (TESTING['MYSTRUCTURE']['MYFARRAYD'], np.float32),
                                                (TESTING['MYSTRUCTURE']['MYDARRAY'], np.float64),
                                                (TESTING['MYSTRUCTURE']['MYDARRAYD'], np.float64),
                                                (TESTING['MYSTRUCTURE']['NESTEDSTRUCT']
                                                 ['MYLARRAY'], np.int32),
                                                (TESTING['MYSTRUCTURE']['NESTEDSTRUCT']
                                                 ['MYLARRAYD'], np.int32),
                                                (TESTING['MYSTRUCTURE']['RANDOMNUMBERS'], np.int16),
                                                (TESTING['MYSTRUCTURE']['MYCARRAY'], complex),
                                                (TESTING['MYSTRUCTURE']['MYDCARRAY'], np.complex64)])
def test_type(myarray, dtype):
    assert myarray.dtype == dtype


def test_date():
    creation_str = TESTING['HEADER']['CREATION']
    creation = datetime.datetime.strptime(creation_str, '%a %b %d %H:%M:%S %Y')
    assert int(''.join(chr(x)
                       for x in TESTING['MYSTRUCTURE']['RANDOMNUMBERS'][-4:])) == creation.year


def test_unpacker():
    """
    This test exercises the Unpacker class that was extracted from the stdlib
    """
    data = b'\x00\x00\x00*\xff\xff\xff\xef\x00\x00\x00\t\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00-?\xf333?\xfeffffff\x00\x00\x00\x0bhello world\x00\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x01\x00\x00\x00\x01\x00\x00\x00\x02\x00\x00\x00\x01\x00\x00\x00\x03\x00\x00\x00\x01\x00\x00\x00\x04\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x04what\x00\x00\x00\x02is\x00\x00\x00\x00\x00\x06hapnin\x00\x00\x00\x00\x00\x06doctor\x00\x00'

    s = b'hello world'
    a = [b'what', b'is', b'hapnin', b'doctor']

    up = genx.Unpacker(data)

    assert up.get_position() == 0

    assert up.unpack_int() == 42
    assert up.unpack_int() == -17
    assert up.unpack_uint() == 9
    assert up.unpack_bool() is True

    # remember position
    pos = up.get_position()
    assert up.unpack_bool() is False

    # rewind and unpack again
    up.set_position(pos)
    assert up.unpack_bool() is False

    assert up.unpack_uhyper() == 45
    np.testing.assert_allclose(up.unpack_float(), 1.9)
    np.testing.assert_allclose(up.unpack_double(), 1.9)
    assert up.unpack_string() == s
    assert up.unpack_list(up.unpack_uint) == list(range(5))
    assert up.unpack_array(up.unpack_string) == a

    up.done()

    with pytest.raises(EOFError):
        up.unpack_uint()
