"""XML helper functions"""

from __future__ import absolute_import
from xml.dom.minidom import parseString #pylint: disable=E0611,F0401

__all__ = ['NotTextNodeError', 'xml_to_dict', 'node_to_dict', 'get_node_text']

#
# Converting XML to a Dictionary
# Author: Christoph Dietze
# URL   : http://code.activestate.com/recipes/116539/
#
class NotTextNodeError(Exception):
    pass

def xml_to_dict(xmlstring):
    """
    Converts an XML string to a Python dictionary

    .. Warning::
        This method does not support multiple inner nodes of the same name but
        with different values.  It always takes the last value.

    Examples
    --------
    ::

        <outer>
            <inner>one</inner>
            <inner>two</inner>
        </outer>

    gives you the dict::

       {u'outer': {u'inner': u'two'}}
    """
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
