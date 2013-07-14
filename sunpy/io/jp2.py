"""JPEG 2000 File Reader"""
from __future__ import absolute_import

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import os
import subprocess
import tempfile
import xml.etree.cElementTree as ET

from matplotlib.image import imread

from sunpy.util.xml import xml_to_dict
from sunpy.map.header import MapHeader

from glymur import Jp2k

__all__ = ['read', 'get_header', 'get_data', 'read_xmlbox', 'which', 'is_float']

def read(filepath, j2k_to_image='opj_decompress'):
    """Reads in the file at the specified location"""
    header = get_header(filepath)
    data = get_data(filepath, j2k_to_image=j2k_to_image)
    
    return data, header 

def get_header(filepath):
    """Reads the header in and saves it as a dictionary"""
    xmlstring = read_xmlbox(filepath, "fits")
    pydict = xml_to_dict(xmlstring)["fits"]
    
    #Fix types
    for k, v in pydict.items():
        if v.isdigit():
            pydict[k] = int(v)
        elif is_float(v):
            pydict[k] = float(v)
            
    # Remove newlines from comment
    if 'comment' in pydict:
        pydict['comment'] = pydict['comment'].replace("\n", "")
            
    return MapHeader(pydict)

def get_data(filepath, j2k_to_image="opj_decompress"):
    """Extracts the data portion of a JPEG 2000 image
    """
    jp2 = Jp2k(filepath)
    data = jp2.read()
    return data

def read_xmlbox(filepath, root):
    """
    Extracts the XML box from a JPEG 2000 image.
    
    Given a filename and the name of the root node, extracts the XML header box
    from a JP2 image.
    """
    jp2 = Jp2k(filepath)
    # Assumes just a single XML box.
    xmlbox = [box for box in jp2.box if box.box_id == 'xml '][0]
    xmlstr = ET.tostring(xmlbox.xml.find('fits'))

    # Fix any malformed XML (e.g. in older AIA data)
    return xmlstr.replace("&", "&amp;")

def is_float(s):
    """Check to see if a string value is a valid float"""
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def which(program):
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

class MissingOpenJPEGBinaryError(OSError):
    """Unable to find OpenJPEG. Please ensure that OpenJPEG binaries are installed in a 
       location within your system's search PATH, or specify the location manually.
    """
    pass
