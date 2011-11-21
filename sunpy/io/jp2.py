"""JPEG 2000 File Reader"""
from __future__ import absolute_import

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import os
import subprocess
import tempfile
from matplotlib.image import imread
from xml.dom.minidom import parseString

def read(filepath):
    """Reads in the file at the specified location"""
    header = get_header(filepath)
    data = get_data(filepath)
    
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
            
    return pydict

def get_data(filepath, j2k_to_image="j2k_to_image"):
    """Extracts the data portion of a JPEG 2000 image
    
    Uses the OpenJPEG j2k_to_image command, if available, to extract the data
    portion of a JPEG 2000 image. The image is first converted to a temporary
    intermediate file (PGM) and then read back in and stored an as ndarray.
    
    NOTE: PIL is also required for Matplotlib to read in PGM images.
    """
    if j2k_to_image == "j2k_to_image" and os.name is "nt":
        j2k_to_image = "j2k_to_image.exe"

    if which(j2k_to_image) is None:
        raise MissingOpenJPEGBinaryError()
    
    jp2filename = os.path.basename(filepath)
    
    tmpname = "".join(os.path.splitext(jp2filename)[0:-1]) + ".pgm"
    tmpfile = os.path.join(tempfile.mkdtemp(), tmpname)
    
    fnull = open(os.devnull, 'w') 
    subprocess.call([j2k_to_image, "-i", filepath, "-o", tmpfile], 
                    stdout=fnull, stderr=fnull)
    fnull.close()
    
    data = imread(tmpfile)    
    os.remove(tmpfile)
    
    return data

def read_xmlbox(filepath, root):
    """
    Extracts the XML box from a JPEG 2000 image.
    
    Given a filename and the name of the root node, extracts the XML header box
    from a JP2 image.
    """
    fp = open(filepath, 'rb')

    xmlstr = ""
    for line in fp:
        xmlstr += line
        if line.find("</%s>" % root) != -1:
            break
    xmlstr = xmlstr[xmlstr.find("<%s>" % root):]
    
    fp.close()

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

#
# Converting XML to a Dictionary
# Author: Christoph Dietze
# URL   : http://code.activestate.com/recipes/116539/
#
class NotTextNodeError(Exception):
    pass

def xml_to_dict(xmlstring):
    """Converts an XML string to a Python dictionary"""
    return node_to_dict(parseString(xmlstring))

def node_to_dict(node):
    """
    node_to_dict() scans through the children of node and makes a dictionary 
    from the content.
    
    Three cases are differentiated:
    1. If the node contains no other nodes, it is a text-node and 
       {nodeName: text} is merged into the dictionary.
    2. If the node has the attribute "method" set to "true", then it's children 
       will be appended to a list and this list is merged to the dictionary in 
       the form: {nodeName:list}.
    3. Else, node_to_dict() will call itself recursively on the nodes children 
       (merging {nodeName: node_to_dict()} to the dictionary).
    """
    dic = {} 
    for n in node.childNodes:
        if n.nodeType != n.ELEMENT_NODE:
            continue
        if n.getAttribute("multiple") == "true":
            # node with multiple children: put them in a list
            l = []
            for c in n.childNodes:
                if c.nodeType != n.ELEMENT_NODE:
                    continue
                l.append(node_to_dict(c))
                dic.update({n.nodeName: l})
            continue
            
        try:
            text = get_node_text(n)
        except NotTextNodeError:
            # 'normal' node
            dic.update({n.nodeName: node_to_dict(n)})
            continue
    
        # text node
        dic.update({n.nodeName: text})
        continue
    return dic

def get_node_text(node):
    """
    scans through all children of node and gathers the text. if node has 
    non-text child-nodes, then NotTextNodeError is raised.
    """
    t = ""
    for n in node.childNodes:
        if n.nodeType == n.TEXT_NODE:
            t += n.nodeValue
        else:
            raise NotTextNodeError
    return t

class MissingOpenJPEGBinaryError(OSError):
    """Unable to find OpenJPEG. Please ensure that OpenJPEG binaries are installed in a 
       location within your system's search PATH, or specify the location manually.
    """
    pass
