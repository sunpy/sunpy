import io
import re

import pytest

from ..mocks import MockObject, MockOpenTextFile, MockHTTPResponse


@pytest.fixture
def mocked_mockobject():
    return MockObject(records=12)


def test_MockObject_illegal_kwargs(mocked_mockobject):
    """
    Any attempt to use a kwarg which has the same name as an attribute/method of the
    underlying object or datastore will raise a ValueError.
    """
    with pytest.raises(ValueError):
        MockObject(records=[], values=1)

    with pytest.raises(ValueError):
        MockObject(items=('a', 'b', 'c'))

    with pytest.raises(ValueError):
        MockObject(__hash__=0x23424)

    # adding a new 'prohibited' attribute will be prevented
    with pytest.raises(ValueError):
        mocked_mockobject['keys'] = [3, 4]


def test_MockObject_attr(mocked_mockobject):
    """
    builtin hasattr & getattr functions, these don't work on dictionaries
    but they do on classes.
    """
    assert hasattr(mocked_mockobject, 'records') is True
    assert hasattr(mocked_mockobject, 'cost') is False

    assert getattr(mocked_mockobject, 'records') == 12

    with pytest.raises(AttributeError):
        getattr(mocked_mockobject, 'jobs')


def test_MockObject_get(mocked_mockobject):
    """
    Getting attributes from `MockObject` using dot and bracket notation.
    """
    assert mocked_mockobject['records'] == 12
    assert mocked_mockobject.records == 12

    with pytest.raises(AttributeError):
        mocked_mockobject.no_key

    with pytest.raises(KeyError):
        mocked_mockobject['not-here']


def test_MockObject_set_get(mocked_mockobject):
    """
    Setting attributes in `MockObject` using bracket notation *not*
    dot notation.
    """

    # Only change the value of existing & new items using 'bracket' notation
    mocked_mockobject['records'] = 45
    assert mocked_mockobject.records == 45
    assert mocked_mockobject['records'] == 45

    # Using 'dot' notation will set a new attribute on 'MockObject' not on the datastore
    # DO NOT DO THIS!
    mocked_mockobject.records = -344

    # This is equivalent to seattr(mocked_mockobject, 'records', -344). Again, don't do this!
    assert mocked_mockobject.records == -344

    # The 'real' value remains unchanged.
    assert mocked_mockobject['records'] == 45


def test_MockObject_len():
    """
    Testing `MockObject.__len__`
    """
    assert len(MockObject(responses=['a', 'b', 'c', 'd'], requests=(1, 2, 3))) == 2


def test_MockObject_del(mocked_mockobject):
    """
    Ensure `MockObject.__delitem__` is *not* implemented.
    """
    with pytest.raises(NotImplementedError):
        del mocked_mockobject['records']


def test_MockObject_iter(mocked_mockobject):
    """
    Test `MockObject.__iter__`
    """
    assert list(iter(mocked_mockobject)) == ['records']


def test_repr_MockObject():
    """
    Test `MockObject.__repr__`
    """
    empty = MockObject()

    mo_p = re.compile(r"^(?P<_><)sunpy\.tests\.mocks\.MockObject \{\} "
                      "at 0x[0-9A-Fa-f]+L?(?(_)>|)$")
    assert mo_p.match(repr(empty)) is not None


def test_read_only_mode_MockOpenTextFile():
    """
    Reading from a read only file, writing should be prohibited.
    """
    new_line = '\n'
    content = r'a{0}bc{0}nd{0}{0}'.format(new_line)

    read_only = MockOpenTextFile('rom.txt', data=content)
    assert read_only.readable() is True
    assert read_only.writable() is False

    with pytest.raises(io.UnsupportedOperation):
        read_only.write('')

    assert read_only.read() == content
    assert read_only.readlines() == ['{0}{1}'.format(line, new_line)
                                     for line in content.split(new_line)]
    read_only.close()

    with pytest.raises(ValueError):
        read_only.readable()

    with pytest.raises(ValueError):
        read_only.writable()

    with pytest.raises(ValueError):
        read_only.read()

    with pytest.raises(ValueError):
        read_only.readlines()


def test_write_only_mode_MockOpenTextFile():
    """
    Writing to to write-only file, reading should be prohibited.
    """
    write_only = MockOpenTextFile('write.txt', 'w')

    assert write_only.readable() is False
    assert write_only.writable() is True

    with pytest.raises(io.UnsupportedOperation):
        write_only.read()

    data = '0123456789'

    num_chars = write_only.write(data)
    assert num_chars == len(data)


def test_read_and_write_MockOpenTextFile():
    """
    Reading & wrtiting to a file with read/write access.
    """
    rd_wr = MockOpenTextFile(mode='r+')

    assert rd_wr.name == 'N/A'
    assert rd_wr.readable() is True
    assert rd_wr.writable() is True

    # Initailly empty
    assert rd_wr.read() == ''

    data = '0123456789'

    num_chars = rd_wr.write(data)
    assert num_chars == len(data)

    assert rd_wr.read() == data

    rd_wr.close()


def test_repr_MockOpenTextFile():
    """
    Test `MockOpenTextFile.__repr__`
    """
    mo_p = re.compile((r"^(?P<_><)sunpy\.tests\.mocks\.MockOpenTextFile file \'a\' "
                       "mode \'r\' at 0x[0-9A-Fa-f]+L?(?(_)>|)$"))

    assert mo_p.match(repr(MockOpenTextFile('a', 'r'))) is not None


def test_MockHTTPResponse():
    """
    Simple tests querying the headers attribute.
    """
    headers = {'Content-Type': 'text/html',
               'Content-Disposition': 'attachment; filename="filename.jpg"'}

    response = MockHTTPResponse(url='http://abc.com', headers=headers)

    assert response.url == 'http://abc.com'

    assert response.headers.get('Content-Disposition') == 'attachment; filename="filename.jpg"'
    assert response.headers.get('Content-Length') is None

    # Key *not* case insensitive
    assert response.headers.get('content-type') is None
