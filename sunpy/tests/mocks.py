"""
Support for creating our own mocked objects

Author: Michael Charlton <m.charlton@mac.com>
"""

import io

from collections import MutableMapping, defaultdict


class MockObject(MutableMapping):
    """
    Object from which we can construct other 'mocked' objects. See
    following examples.

    Limitations
    -----------

    On initiation a `ValueError` will be raised if any of the `kwargs` have the same name as an
    existing attribute/method of the `MockObject` or underlying data store.

    Updating existing attributes, or adding new ones, should only be done using bracket
    notation and *not* dot notation. Using dot notation will update the `MockObject` and not the
    data store.

    Examples
    --------
    >>> mo = MockObject(code=400)
    >>> mo.code
    400
    >>> mo['code']
    400
    >>> mo.code = 210
    >>> mo.code
    210
    >>> mo['code']
    400

    The recommended way of changing the value of an existing, or new, attribute is using bracket
    notation:

    >>> m = MockObject(start='now')
    >>> m['start'] = 'Thursday'
    >>> m['start']
    'Thursday'
    >>> m.start
    'Thursday'
    """
    def __init__(self, *args, **kwargs):
        self._datastore = dict()
        self.prohibited_attrs = set(dir(self))
        self.prohibited_attrs.update(dir(self._datastore))

        for candidate in kwargs.keys():
            if candidate in self.prohibited_attrs:
                raise ValueError("kwarg '{kwarg}' is already an attribute "
                                 "of {datastore} or {obj}".format(kwarg=candidate,
                                                                  datastore=type(self._datastore),
                                                                  obj=type(self)))
        self.update(dict(*args, **kwargs))

    def __getattr__(self, name):
        if name in self._datastore:
            return self._datastore[name]

        raise AttributeError(name)

    def __getitem__(self, name):
        return self._datastore[name]

    def __setitem__(self, name, value):
        if name in self.prohibited_attrs:
            raise ValueError("Name '{name}' is already an attribute "
                             "of {datastore} or {obj}".format(name=name,
                                                              datastore=type(self._datastore),
                                                              obj=type(self)))
        self._datastore[name] = value

    def __delitem__(self, name):
        raise NotImplementedError("'del' operation for {0} "
                                  "not supported".format(self.__class__.__name__))

    def __iter__(self):
        return iter(self._datastore)

    def __len__(self):
        return len(self._datastore)

    def __repr__(self):
        return ("<{module}.{name} {contents} at {address}>".format(module=self.__module__,
                                                                   name=self.__class__.__name__,
                                                                   contents=self._datastore,
                                                                   address=hex(id(self))))


class MockHTTPResponse(MockObject):
    """
    The result of calling `urlopen(...)`. For this implementation we are only interested in
    querying the 'headers' attribute, which is a http.client.HTTPMessage object.

    Parameters
    ----------

    url : `str` optional, default: ''
        The url of the connection

    headers : `dict` of `str` optional, default: empty dictionary
        HTTP header fields of the response message.

    Limitations
    -----------

    On a 'real' http.client.HTTPMessage, header name retrieval is case insensitive.
    In this implementation the header names are case sensitive.

    Examples
    --------
    >>> result = MockHTTPResponse(url='http://abc.com', headers={'Content-Type':'text/html'})
    >>> result.headers.get('Content-Type')
    'text/html'
    """
    def __init__(self, *args, **kwargs):
        super(MockHTTPResponse, self).__init__(*args, **kwargs)

        self.setdefault('url', '')

        headers_store = defaultdict(lambda: None)

        if 'headers' in self:
            headers_store.update(self['headers'])

        self['headers'] = headers_store


class MockOpenTextFile(MockObject):
    """
    Partial implementation of a file like object for reading/wrtiing text files. Binary files
    are *not* supported.

    Many methods not implemented, no attempt is made to keep track of where we are
    in the file say in regards to any read operation.

    Parameters
    ----------
    file : `str` optional, default:'UNKNOWN'
        The name of the file. As per the builtin `open` function

    mode : `str` optional, default:'r'
        The way in which the file is to be used. As per the builtin `open` function

    data : `str` optional keyword, default:''
        The inital data which can be read from the file.

    Limitations
    -----------

    Unlike in a real file, this implementation makes no attempt to keep
    track of where we are, when reading or writing.

    Examples
    --------

    >>> dummy_read_only = MockOpenTextFile()
    >>> named_write = MockOpenTextFile(file='a.txt', mode='w')
    >>> named_read = MockOpenTextFile('b.txt')
    >>> named_rd_wr = MockOpenTextFile('c.txt', 'r+', data='Hello, world')
    """
    def __init__(self, *args, **kwargs):

        # Positional and/or keword args can be used for the 'file' & 'mode'
        # parameters. Could do a lot more checking to make sure all required
        # arguments are present

        num_pos_args = len(args)

        if num_pos_args == 1:
            kwargs['file'] = args[0]
        elif num_pos_args == 2:
            kwargs['file'] = args[0]
            kwargs['mode'] = args[1]

        if 'mode' not in kwargs:
            kwargs['mode'] = 'r'

        super(MockOpenTextFile, self).__init__(**kwargs)
        self.setdefault('file', 'N/A')
        self['name'] = self['file']
        self.setdefault('closed', False)
        self.setdefault('data', '')

    def write(self, content):
        if not self.writable():
            raise io.UnsupportedOperation(':not writable')

        self.data += content

        return len(content)

    def read(self):
        if not self.readable():
            raise io.UnsupportedOperation(': not readable')

        return self.data

    def readlines(self):
        if self.closed:
            raise ValueError('I/O operation on closed file')

        # Documentation recommends using '\n' as the line terminator when reading/writing text
        # files. See `os.linesep` in https://docs.python.org/3/library/os.html
        new_line = '\n'
        return ['{0}{1}'.format(line, new_line) for line in self.data.split(new_line)]

    def readable(self):
        if self.closed:
            raise ValueError('I/O operation on closed file')

        return 'r' in self.mode

    def writable(self):
        if self.closed:
            raise ValueError('I/O operation on closed file')

        return ('w' in self.mode) or ('r+' in self.mode)

    def close(self):
        self.closed = True
        self.data = ''

    def __repr__(self):
        return ("<{module}.{name} file '{file}' mode '{mode}' "
                "at {address}>".format(module=self.__module__, name=self.__class__.__name__,
                                       file=self.file, mode=self.mode,
                                       address=hex(id(self))))
