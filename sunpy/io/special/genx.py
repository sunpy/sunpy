"""
This module implements a solarsoft genx file reader.
"""
import copy
import struct
from collections import OrderedDict

import numpy as np

__all__ = ['read_genx']


# This class has been copied from the Python 3.11 stdlib xdrlib.py file under
# the terms of the PSF licence 2.0
class Error(Exception):
    """Exception class for this module. Use:
    except xdrlib.Error as var:
        # var has the Error instance for the exception
    Public ivars:
        msg -- contains the message
    """

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        return repr(self.msg)

    def __str__(self):
        return str(self.msg)


class ConversionError(Error):
    pass


class Unpacker:
    """Unpacks various data representations from the given buffer."""

    def __init__(self, data):
        self.reset(data)

    def reset(self, data):
        self.__buf = data
        self.__pos = 0

    def get_position(self):
        return self.__pos

    def set_position(self, position):
        self.__pos = position

    def get_buffer(self):
        return self.__buf

    def done(self):
        if self.__pos < len(self.__buf):
            raise Error('unextracted data remains')

    def unpack_uint(self):
        i = self.__pos
        self.__pos = j = i+4
        data = self.__buf[i:j]
        if len(data) < 4:
            raise EOFError
        return struct.unpack('>L', data)[0]

    def unpack_int(self):
        i = self.__pos
        self.__pos = j = i+4
        data = self.__buf[i:j]
        if len(data) < 4:
            raise EOFError
        return struct.unpack('>l', data)[0]

    unpack_enum = unpack_int

    def unpack_bool(self):
        return bool(self.unpack_int())

    def unpack_uhyper(self):
        hi = self.unpack_uint()
        lo = self.unpack_uint()
        return int(hi) << 32 | lo

    def unpack_hyper(self):
        x = self.unpack_uhyper()
        if x >= 0x8000000000000000:
            x = x - 0x10000000000000000
        return x

    def unpack_float(self):
        i = self.__pos
        self.__pos = j = i+4
        data = self.__buf[i:j]
        if len(data) < 4:
            raise EOFError
        return struct.unpack('>f', data)[0]

    def unpack_double(self):
        i = self.__pos
        self.__pos = j = i+8
        data = self.__buf[i:j]
        if len(data) < 8:
            raise EOFError
        return struct.unpack('>d', data)[0]

    def unpack_fstring(self, n):
        if n < 0:
            raise ValueError('fstring size must be nonnegative')
        i = self.__pos
        j = i + (n+3)//4*4
        if j > len(self.__buf):
            raise EOFError
        self.__pos = j
        return self.__buf[i:i+n]

    unpack_fopaque = unpack_fstring

    def unpack_string(self):
        n = self.unpack_uint()
        return self.unpack_fstring(n)

    unpack_opaque = unpack_string
    unpack_bytes = unpack_string

    def unpack_list(self, unpack_item):
        alist = []
        while (x := self.unpack_uint()) != 0:
            if x != 1:
                raise ConversionError(f'0 or 1 expected, got {x!r}')
            item = unpack_item()
            alist.append(item)
        return alist

    def unpack_farray(self, n, unpack_item):
        return [unpack_item() for _ in range(n)]

    def unpack_array(self, unpack_item):
        n = self.unpack_uint()
        return self.unpack_farray(n, unpack_item)


class SSWUnpacker(Unpacker):
    """
    A `xdrlib.Unpacker` customisation to read strings and complex data as
    written by IDL.
    """

    def unpack_string(self):
        n = self.unpack_uint()
        if n > 0:
            n = self.unpack_uint()
        return self.unpack_fstring(n).decode('utf-8')

    def unpack_complex(self):
        return complex(self.unpack_float(), self.unpack_float())

    def unpack_complex_double(self):
        return complex(self.unpack_double(), self.unpack_double())


def read_struct_skeleton(xdrdata):
    """
    Reads the skeleton of the IDL's structure as written in solarsoft
    ``build_str()`` function.
    """
    # Read Number of tags
    ntags = xdrdata.unpack_uint()
    # Read names of tags
    tags = [xdrdata.unpack_string() for n in range(ntags)]
    # Read the tag type
    tagdict = OrderedDict()
    for tt in tags:
        dim = xdrdata.unpack_uint()
        arr_size = xdrdata.unpack_farray(dim + 2, xdrdata.unpack_int)  # [7, 1]
        typedata = arr_size[-2]
        if typedata == 8:  # it's a structure
            if dim == 2 and arr_size[0] == 1:
                # For when structures has been defined with 2 dim but only has one:
                # bb = replicate({tata: 1, bebe:2, casa:'asa'}, 3)
                # dd = replicate({tata: 1, bebe:2, casa:'asa'}, 1, 3)
                # print,size(bb), size(dd)
                #            1           3           8           3
                #            2           1           3           8           3
                dim = 1
                arr_size[0] = arr_size[1]
            if arr_size[-1] > 1:
                tagdict[tt] = np.array([read_struct_skeleton(xdrdata)] *
                                       arr_size[-1]).reshape(arr_size[0:-2])
            else:
                tagdict[tt] = read_struct_skeleton(xdrdata)
        else:
            tagdict[tt] = [dim] + arr_size
    return tagdict


def struct_to_data(xdrdata, subskeleton):
    """
    Converts the dictionary with the keys and IDL's size output to the data
    stored in the xdrdata.

    ``subskeleton`` must contain the size and type of the data that's going to be read in the right
    order (that's why `OrderedDict` is used). Then the data is read and the ``subskeleton`` is updated
    with the data itself.

    Parameters
    ----------
    xdrdata : The return from `~sunpy.io.special.SSWUnpacker`
        The data to read in from `~sunpy.io.special.SSWUnpacker`.
    subskeleton : `OrderedDict` or `np.array`
        Contains the size and type of the ``xdrdata``.
    """
    # http://www.harrisgeospatial.com/docs/SIZE.html
    types_dict = {
        2: (xdrdata.unpack_int, np.int16),  # int
        3: (xdrdata.unpack_int, np.int32),  # long
        4: (xdrdata.unpack_float, np.float32),
        5: (xdrdata.unpack_double, np.float64),
        6: (xdrdata.unpack_complex, complex),
        7: (xdrdata.unpack_string, None),
        9: (xdrdata.unpack_complex_double, np.complex64),
        12: (xdrdata.unpack_uint, np.uint16),  # unsign int
        13: (xdrdata.unpack_uint, np.uint32),  # unsign Long int
        14: (xdrdata.unpack_hyper, np.int64),  # Long64
        15: (xdrdata.unpack_uhyper, np.uint64),  # unsign Long64
    }
    for key in subskeleton:
        if isinstance(subskeleton[key], OrderedDict):
            struct_to_data(xdrdata, subskeleton[key])
        elif isinstance(subskeleton[key], np.ndarray):
            testlist = list()
            struct_shape = subskeleton[key].shape
            for elem in subskeleton[key].flatten():
                elem2 = copy.deepcopy(elem)
                struct_to_data(xdrdata, elem2)
                testlist.append(elem2)
            subskeleton[key] = np.array(testlist).reshape(struct_shape)
        else:
            sswsize = subskeleton[key]
            sswtype = sswsize[-2]
            if sswsize[0] == 0:
                subskeleton[key] = types_dict[sswtype][0]()
            else:
                subskeleton[key] = np.array(xdrdata.unpack_farray(sswsize[-1], types_dict[sswtype][0]),
                                            dtype=types_dict[sswtype][1]).reshape(sswsize[1:-2][::-1])


def read_genx(filename):
    """
    solarsoft genx file reader.

    genx files have been used to store calibration data for multiple
    instruments and distributed within solarsoft. They are stored in XDR
    format; The External Data Representation Standard file format (XDR) is
    described in `RFC 1014 <https://www.rfc-editor.org/rfc/rfc1014.txt>`__,
    written by Sun Microsystems, Inc. June 1987.

    SolarSoft genx writer creates structures to store the values together with
    the variable names. It use the ``size`` IDL function to include the data
    type, dimension and number of elements that each variable contains.

    Parameters
    ----------
    filename : `str`
        The genx file to be read

    Returns
    -------
    output : `~collections.OrderedDict`
        A dictionary with possibly nested dictionaries with the data in
        the genx file.

    Notes
    -----
    The reader aims to maintain the shape and type of the arrays, but take care with the
    difference in indexing between Python and IDL (row mayor vs column mayor).

    Regarding the type notice that single numbers are converted to python precision, therefore
    a single integer is converted from 16 to 32/64 bits, and a float from 32 to 64.

    **Strings** read from genx files are assumed to be UTF-8.
    """
    with open(filename, mode='rb') as xdrfile:
        xdrdata = SSWUnpacker(xdrfile.read())

    # HEADER information
    version, xdr = xdrdata.unpack_int(), xdrdata.unpack_int()
    creation = xdrdata.unpack_string()
    if version == 2:
        arch = xdrdata.unpack_string()
        os = xdrdata.unpack_string()
        release = xdrdata.unpack_string()
    text = xdrdata.unpack_string()

    # TODO: I don't think will have dim>1 but need
    # to check, if it's larger like savegen has run
    # with a multidimensional structure, then
    # probably the read_struct_skeleton will have to
    # run as in multi-dim structure.

    dim = xdrdata.unpack_int()
    # [1, 8, 1] = Main structure for the data
    arr_size = xdrdata.unpack_farray(dim + 2, xdrdata.unpack_int)
    # the number of upper level strs
    # This is used somehow
    mainsize = arr_size[2]  # NOQA

    skeleton = read_struct_skeleton(xdrdata)
    struct_to_data(xdrdata, skeleton)
    xdrdata.done()
    skeleton['HEADER'] = OrderedDict([('VERSION', version), ('XDR', xdr), ('CREATION', creation)])
    if version == 2:
        skeleton['HEADER']['IDL_VERSION'] = OrderedDict([('ARCH', arch),
                                                         ('OS', os),
                                                         ('RELEASE', release)])
    skeleton['HEADER']['TEXT'] = text
    skeleton.move_to_end('HEADER', last=False)
    return skeleton
