import datetime

import numpy as np
import pytest

from sunpy.data.test import get_test_filepath
from sunpy.io.special import genx

TESTING = genx.read_genx(get_test_filepath('generated_sample.genx'))
TESTING_BYTE = genx.read_genx(get_test_filepath('byte_sample.genx'))
TESTING_BYTE_PADDING = genx.read_genx(get_test_filepath('byte_padding_sample.genx'))


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


def test_byte_data():
    # Byte data is written by IDL as a single counted block (the xdr_bytes
    # representation), not element by element like the other numeric types.
    # The sample file holds a byte array, a scalar byte, a 2D byte array
    # and an int array stored after them.
    np.testing.assert_array_equal(TESTING_BYTE['MYBYTE_ARRAY'],
                                  np.array([0, 1, 2, 127, 128, 255], dtype=np.uint8))
    assert TESTING_BYTE['MYBYTE_ARRAY'].dtype == np.uint8
    assert TESTING_BYTE['MYBYTE'] == 200
    np.testing.assert_array_equal(TESTING_BYTE['MYBYTE_ARRAY_DIMENSION'],
                                  np.arange(6, dtype=np.uint8).reshape(3, 2))
    # Reading the byte blocks must consume exactly the right number of
    # bytes, so the int array stored after them is not corrupted.
    np.testing.assert_array_equal(TESTING_BYTE['MYINT_ARRAY'],
                                  np.array([0, 1, 2], dtype=np.int16))
    # Byte data nested inside a structure.
    np.testing.assert_array_equal(TESTING_BYTE['MYSTRUCTURE']['MYBYTES'],
                                  np.array([10, 20, 30, 40], dtype=np.uint8))
    assert TESTING_BYTE['MYSTRUCTURE']['MYINT'] == 7


@pytest.mark.parametrize(('tag', 'values', 'following'), [
    ('PAD0', [], 101),
    ('PAD1', [11], 102),
    ('PAD2', [21, 22], 103),
    ('PAD3', [31, 32, 33], 104),
    ('PAD4', [41, 42, 43, 44], 105),
])
def test_byte_array_padding(tag, values, following):
    # A byte array is stored as a counted block padded up to a multiple
    # of 4 bytes, so lengths 0 to 4 exercise every padding width. Reading
    # must consume the block and its padding exactly, which is checked
    # against the int value stored after each array.
    array = TESTING_BYTE_PADDING[tag]
    np.testing.assert_array_equal(array, np.array(values, dtype=np.uint8))
    assert array.dtype == np.uint8
    assert TESTING_BYTE_PADDING[f'AFTER{tag[-1]}'] == following


def test_byte_scalar_empty_block_error():
    # A scalar byte whose counted block is empty is malformed.
    with pytest.raises(ValueError, match='scalar byte'):
        genx.read_genx(get_test_filepath('byte_scalar_empty_sample.genx'))


def test_byte_scalar_long_block_error():
    # A scalar byte whose counted block holds more than one byte is
    # malformed; the extra bytes must not be silently discarded.
    with pytest.raises(ValueError, match='scalar byte'):
        genx.read_genx(get_test_filepath('byte_scalar_long_sample.genx'))


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


def test_read_genx_error(tmp_path):
    idl_signature = b'SR\x00\x04' + b'\x00' * 10
    temp_file = tmp_path / "dummy_idl.geny"
    temp_file.write_bytes(idl_signature)
    with pytest.raises(ValueError, match="is an IDL savefile, not a genx file"):
        genx.read_genx(temp_file)


def test_read_genx_too_short_error(tmp_path):
    short_file = b'12345'
    temp_file = tmp_path / "short_dummy.genx"
    temp_file.write_bytes(short_file)
    with pytest.raises(ValueError, match="is too short to be a valid genx file"):
        genx.read_genx(temp_file)


def test_read_genx_unsupported_version_error(tmp_path):
    text_data = b'this is a plain text file'
    temp_file = tmp_path / "plain_text.txt"
    temp_file.write_bytes(text_data)
    with pytest.raises(ValueError, match="is not a 'genx' file with version"):
        genx.read_genx(temp_file)
