import xdrlib
from collections import OrderedDict
from functools import partial
import copy

import numpy as np

__all__ = ['read_genx']

class SSWUnpacker(xdrlib.Unpacker):
    """
    `xdrlib.Unpacker` customisation to read strings and complex data as written
    by IDL.
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
    Reads the skeleton of the IDL's structure as written in
    solarsoft `build_str()` function.
    """
    # Read Number of tags
    ntags = xdrdata.unpack_uint()
    # Read names of tags
    tags = [xdrdata.unpack_string() for n in range(ntags)]
    # Read the tag type
    tagdict = OrderedDict()
    for tt in tags:
        dim = xdrdata.unpack_uint()
        arr_size = xdrdata.unpack_farray(dim + 2, xdrdata.unpack_int) # [7, 1]
        typedata = arr_size[-2]
        if typedata == 8: # it's a structure
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
                tagdict[tt] = np.array([read_struct_skeleton(xdrdata)] * arr_size[-1]).reshape(arr_size[0:-2])
            else:
                tagdict[tt] = read_struct_skeleton(xdrdata)
        else:
            tagdict[tt] = [dim] + arr_size
    return tagdict

def struct_to_data(xdrdata, subskeleton):
    """"
    Converts the dictionary with the keys and IDL's size output to
    the data stored in the xdrdata.

    `subskeleton` must contain the size and type of the data that's going to be
    read in the right order (that's why `OrderedDict` is used). Then the data is
    read and the `subskeleton` is updated with the data itself.
    """
    #http://www.harrisgeospatial.com/docs/SIZE.html
    types_dict = {
        2: (xdrdata.unpack_int, np.int16), # int
        3: (xdrdata.unpack_int, np.int32), # long
        4: (xdrdata.unpack_float, np.float32),
        5: (xdrdata.unpack_double, np.float64),
        6: (xdrdata.unpack_complex, np.complex),
        7: (xdrdata.unpack_string, None),
        9: (xdrdata.unpack_complex_double, np.complex64),
        12: (xdrdata.unpack_uint, np.uint16), # unsign int
        13: (xdrdata.unpack_uint, np.uint32), # unsign Long int
        14: (xdrdata.unpack_hyper, np.int64), # Long64
        15: (xdrdata.unpack_uhyper, np.uint64), # unsign Long64
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
    """solarsoft genx file reader

    genx files have been used to store calibration data for multiple
    instruments and distributed within solarsoft. They are stored in XDR
    format; The External Data Representation Standard file format (XDR) is
    described in `RFC 1014 <https://tools.ietf.org/html/rfc1014>`_,
    written by Sun Microsystems, Inc. June 1987

    SolarSoft genx writer creates structures to store the values together with
    the variable names. It use the `size` IDL function to include the data
    type, dimension and number of elements that each variable contains.

    Parameters
    ----------
    filename : `str`
        The genx file to be read

    Returns
    -------
    output : `OrderedDict`
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
    creation =  xdrdata.unpack_string()
    if version == 2:
        arch =  xdrdata.unpack_string()
        os =  xdrdata.unpack_string()
        release =  xdrdata.unpack_string()
    text = xdrdata.unpack_string()


    dim = xdrdata.unpack_int() # TODO: I don't think will have dim>1 but need
                               # to check, if it's larger like savegen has run
                               # with a multidimensional structure, then
                               # probably the read_struct_skeleton will have to
                               # run as in multi-dim structure.
    arr_size = xdrdata.unpack_farray(dim + 2, xdrdata.unpack_int) # [1, 8, 1] = Main structure for the data
    mainsize = arr_size[2] # number of upper level strs
    skeleton = read_struct_skeleton(xdrdata)
    struct_to_data(xdrdata, skeleton)
    xdrdata.done()
    skeleton['HEADER'] = OrderedDict([('VERSION', version), ('XDR', xdr), ('CREATION', creation)])
    if version == 2:
        skeleton['HEADER']['IDL_VERSION'] = OrderedDict([('ARCH', arch),
                                                         ('OS', os),
                                                         ('RELEASE', release)])
    skeleton['HEADER']['TEXT'] = text
    # TODO: for python >= 3.2; so we can keep the original order as how it's stored in the file
    # skeleton.move_to_end('HEADER', last=False)
    return skeleton
