from __future__ import absolute_import
"""
JPEG 2000 File Reader
"""
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
    xmldict = xml_to_dict(xmlstring)["fits"]
    
    #Fix types
    for k, v in xmldict.items():
        if v.isdigit():
            xmldict[k] = int(v)
        elif is_float(v):
            xmldict[k] = float(v)
            
    return xmldict

def get_data(filepath):
    """Extracts the data portion of a JPEG 2000 image
    
    Uses the OpenJPEG j2k_to_image command, if available, to extract the data
    portion of a JPEG 2000 image. The image is first converted to a temporary
    intermediate file (PGM) and then read back in and stored an as ndarray.
    """
    jp2filename = os.path.basename(filepath)
    
    tmpname = "".join(os.path.splitext(jp2filename)[0:-1]) + ".pgm"
    tmpfile = os.path.join(tempfile.mkdtemp(), tmpname)
    
    fnull = open(os.devnull, 'w') 
    subprocess.call(["j2k_to_image", "-i", filepath, "-o", tmpfile], 
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
