import xdrlib
from collections import OrderedDict
from functools import partial

class ssw_xdrlib(xdrlib.Unpacker):
    """
    `xdrlib.Unpacker` customisation to read strings as written by IDL.
    """
    def unpack_string(self):
        n = self.unpack_uint()
        if n > 0:
            n = self.unpack_uint()
            print('second read', n)
        return self.unpack_fstring(n).decode()

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
        if typedata == 8: # it's not a structure
            if dim == 2 and arr_size[0] == 1:
                # For when structures has been defined with 2 dim but only has one:
                # bb = replicate({tata: 1, bebe:2, casa:'asa'}, 3)
                # dd = replicate({tata: 1, bebe:2, casa:'asa'}, 1, 3)
                # print,size(bb), size(dd)
                #            1           3           8           3
                #            2           1           3           8           3
                dim = 1
                arr_size[0] = arr_size[1]
            tagdict[tt] = read_struct_skeleton(xdrdata)
        else:
            tagdict[tt] = [dim] + arr_size
    return tagdict

def unsuported_dtype(dtype):
    raise ValueError("Can't read {} data type".format(dtype))

def struct_to_data(xdrdata, subskeleton):
    """"
    Converts the dictionary with the keys and IDL's size output to
    the data stored in the xdrdata.
    """
    #http://www.harrisgeospatial.com/docs/SIZE.html
    types_dict = {
        2: xdrdata.unpack_int, # int
        3: xdrdata.unpack_int, # long
        4: xdrdata.unpack_float,
        5: xdrdata.unpack_double,
        6: partial(unsuported_dtype, 'Complex'), # FIXME: How to read complex?
        7: xdrdata.unpack_string,
        9: partial(unsuported_dtype, 'Double Complex'), # FIXME: How to read Double complex?
        12: xdrdata.unpack_uint, # unsign int
        13: xdrdata.unpack_uint, # unsign Long int
        14: partial(unsuported_dtype, 'Long 64'), # Long64          # FIXME: How to do 64?
        15: partial(unsuported_dtype, 'Unsign Long 64'), # unsign Long64   # FIXME: How to do 64?
    } 
    for key in subskeleton:
        if isinstance(subskeleton[key], OrderedDict):
            struct_to_data(xdrdata, subskeleton[key])
        else:
            sswsize = subskeleton[key]
            sswtype = sswsize[-2]
            if sswsize[0] == 0:
                subskeleton[key] = types_dict[sswtype]()
            else:
                subskeleton[key] = xdrdata.unpack_farray(sswsize[-1], types_dict[sswtype])

def read_genx(filename):
    """
    solarsoft genx file reader

    Parameters
    ----------
    filename : `str`
        The genx file to be read

    Returns
    -------
    output : `OrderedDict`
        A dictionary with possibly nested dictionaries with the data in
        the genx file.

    """
    # TODO: check filename exists
    with open(filename, mode='rb') as xdrfile:
        xdrdata = ssw_xdrlib(xdrfile.read())

    # HEADER information
    version, xdr = xdrdata.unpack_int(), xdrdata.unpack_int() # TODO: what when version 1? #TODO: what when not xdr?
    creation =  xdrdata.unpack_string()
    arch =  xdrdata.unpack_string()
    os =  xdrdata.unpack_string()
    release =  xdrdata.unpack_string()
    text = xdrdata.unpack_string()


    dim = xdrdata.unpack_int() # so it gives 1 (as in gdl) # TODO: I don't think will have dim>1 but need to check
    arr_size = xdrdata.unpack_farray(dim + 2, xdrdata.unpack_int) # [1, 8, 1] = Main structure for the data
    mainsize = arr_size[2] # number of upper level strs
    skeleton = read_struct_skeleton(xdrdata)
    struct_to_data(xdrdata, skeleton)
    xdrdata.done()
    skeleton['HEADER'] = OrderedDict([('version', version), ('xdr', xdr), ('creation', creation),
                                      ('idl_version', OrderedDict([('arch', arch),
                                                                   ('os', os),
                                                                   ('release', release)])),
                                      ('text', text)])
    # TODO: for python >= 3.2; so we can keep the original order as how it's stored in the file
    # skeleton.move_to_end('HEADER', last=False)
    return skeleton
