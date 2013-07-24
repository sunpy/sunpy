"""JPEG 2000 File Reader"""
from __future__ import absolute_import

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import os
import xml.etree.cElementTree as ET

from glymur import Jp2k
import glymur.jp2box

from sunpy.util.xml import xml_to_dict
from sunpy.io.header import FileHeader

__all__ = ['read', 'get_header', 'write']

def read(filepath):
    """
    Reads a JPEG2000 file
    
    Parameters
    ----------
    filepath : string
        The file to be read
    
    j2k_to_image : string
        binary to use for reading?
    
    Returns
    -------
    pairs : list
        A list of (data, header) tuples
    """
    header = get_header(filepath)[0]
    data = _get_data(filepath)
    
    return [(data, header)]

def get_header(filepath):
    """
    Reads the header form the file
    
    Parameters
    ----------
    filepath : string
        The file to be read
        
    Returns
    -------
    headers : list
        A list of headers read from the file
    """
    xmlstring = _read_xmlbox(filepath, "fits")
    pydict = xml_to_dict(xmlstring)["fits"]
    
    #Fix types
    for k, v in pydict.items():
        if v.isdigit():
            pydict[k] = int(v)
        elif _is_float(v):
            pydict[k] = float(v)
            
    # Remove newlines from comment
    if 'comment' in pydict:
        pydict['comment'] = pydict['comment'].replace("\n", "")
            
    return [FileHeader(pydict)]

def write(fname, data, header):
    """
    Place holder for required file writer
    """
    raise NotImplementedError("No jp2 writer is implemented")

def _get_data(filepath):
    """Extracts the data portion of a JPEG 2000 image
    
    The image as read back in is upside down, and so it is spun around for the
    correct orientation. -- Is this true???
    """
    jp2 = Jp2k(filepath)
    data = jp2.read()
    return data

def _read_xmlbox(filepath, root):
    """
    Extracts the XML box from a JPEG 2000 image.
    
    Given a filename and the name of the root node, extracts the XML header box
    from a JP2 image.
    """
    jp2 = Jp2k(filepath)
    # Assumes just a single XML box.
    xmlbox = [box for box in jp2.box if box.box_id == 'xml '][0]
    xmlstr = ET.tostring(xmlbox.xml.find(root))

    # Fix any malformed XML (e.g. in older AIA data)
    return xmlstr.replace("&", "&amp;")

def _is_float(s):
    """Check to see if a string value is a valid float"""
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def _which(program):
    """Checks for existence of executable
    
    Source: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python/377028#377028
    """
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program) #pylint: disable=W0612

    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
