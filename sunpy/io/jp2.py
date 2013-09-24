"""JPEG 2000 File Reader"""
from __future__ import absolute_import

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import os
import subprocess
import tempfile

from matplotlib.image import imread

from sunpy.util.xml import xml_to_dict
from sunpy.io.header import FileHeader

__all__ = ['read', 'get_header', 'write']

def read(filepath, j2k_to_image='j2k_to_image'):
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
    header = get_header(filepath)
    data = _get_data(filepath, j2k_to_image=j2k_to_image)
    
    return [(data, header[0])]

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

def _get_data(filepath, j2k_to_image="opj_decompress"):
    """Extracts the data portion of a JPEG 2000 image
    
    Uses the OpenJPEG j2k_to_image command, if available, to extract the data
    portion of a JPEG 2000 image. The image is first converted to a temporary
    intermediate file (PNG) and then read back in and stored an as ndarray.
    
    The image as read back in is upside down, and so it is spun around for the
    correct orientation.
    """
    if os.name is "nt":
        if (j2k_to_image == "j2k_to_image") or (j2k_to_image == "opj_decompress"):
            j2k_to_image = j2k_to_image+".exe"

    if _which(j2k_to_image) is None:
        raise MissingOpenJPEGBinaryError("You must first install the OpenJPEG "
                                         "(version >=1.4) binaries before using "
                                         "this functionality.")
    
    jp2filename = os.path.basename(filepath)
    
    tmpname = "".join(os.path.splitext(jp2filename)[0:-1]) + ".png"
    tmpfile = os.path.join(tempfile.mkdtemp(), tmpname)
    
    with open(os.devnull, 'w') as fnull:
        subprocess.call([j2k_to_image, "-i", filepath, "-o", tmpfile], 
                        stdout=fnull, stderr=fnull)
    
    data = imread(tmpfile)    
    os.remove(tmpfile)
    
    # flip the array around since it has been read in upside down.
    return data[::-1]

def _read_xmlbox(filepath, root):
    """
    Extracts the XML box from a JPEG 2000 image.
    
    Given a filename and the name of the root node, extracts the XML header box
    from a JP2 image.
    """
    with open(filepath, 'rb') as fp:

        xmlstr = ""
        for line in fp:
            xmlstr += line
            if line.find("</%s>" % root) != -1:
                break

        start = xmlstr.find("<%s>" % root)
        end = xmlstr.find("</%s>" % root) + len("</%s>" % root)
        
        xmlstr = xmlstr[start : end]

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

class MissingOpenJPEGBinaryError(OSError):
    """Unable to find OpenJPEG. Please ensure that OpenJPEG binaries are installed in a 
       location within your system's search PATH, or specify the location manually.
    """
    pass
