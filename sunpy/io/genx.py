
import xdrlib
from collections import OrderedDict

def read_struct_skeleton(xdrdata):
    # Read Number of tags
    ntags = xdrdata.unpack_uint()
    # Read names of tags
    tags = list()
    for n in range(ntags):
        len_str = xdrdata.unpack_int()
        if len_str > 0:
            tags.append(xdrdata.unpack_string().decode())
        else:
            tags.append('')
    # Read the tag type
    tagdict = OrderedDict()
    for tt in tags:
        dim = xdrdata.unpack_uint()
        arr_size = xdrdata.unpack_farray(dim + 2, xdrdata.unpack_int) # [7, 1]
        typedata = arr_size[-2]
        if typedata == 8: # it's not a structure
            if dim == 2 and arr_size[0] == 1:
                # For when struct has been defined with 2 dim but only has one:
                # bb = replicate({tata: 1, bebe:2, casa:'asa'}, 3)
                dim = 1
                arr_size[0] = arr_size[1]
            tagdict[tt] = read_str_skeleton(xdrdata)
        else:
            tagdict[tt] = [dim] + arr_size
    return tagdict


def read_genx(filename):
    filename = './aia_V6_all_fullinst.genx'
    with open(filename, mode='rb') as xdrfile:
        xdrdata = xdrlib.Unpacker(xdrfile.read())
    version, xdr = xdrdata.unpack_int(), xdrdata.unpack_int()
    creation =  xdrdata.unpack_uint()
    creation =  xdrdata.unpack_string() # they are 24, but padded to 4 byte alignment.
    archN =  xdrdata.unpack_uint()
    arch =  xdrdata.unpack_string()
    osN =  xdrdata.unpack_uint()
    os =  xdrdata.unpack_string()
    releaseN =  xdrdata.unpack_uint()
    release =  xdrdata.unpack_string()
    textn = xdrdata.unpack_uint()
    text = xdrdata.unpack_fstring(textn)
    #text = xdrdata.unpack_fstring(textn) # in this case there's not

    dim = xdrdata.unpack_int() # so it gives 1 (as in gdl)
    arr_size = xdrdata.unpack_farray(dim + 2, xdrdata.unpack_int) # [1, 8, 1]
    mainsize = arr_size[2] # number of upper level strs
    return read_str_skeleton(xdrdata)
