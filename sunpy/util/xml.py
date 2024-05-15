"""
This module provides XML helper functions.
"""
from xml.dom.minidom import parseString

__all__ = ['NotTextNodeError', 'xml_to_dict', 'node_to_dict', 'get_node_text']


class NotTextNodeError(Exception):
    pass


def xml_to_dict(xmlstring):
    """
    Converts an XML string to a Python dictionary.

    .. warning::
        This method does not support multiple inner nodes of the same name but
        with different values. It always takes the last value.

    Parameters
    ----------
    xmlstring : `str`
        A `str` of xml content.

    Returns
    -------
    `dict`
        The string xml input as a dictionary.

    Examples
    --------

    .. code-block:: xml

        <outer>
            <inner>one</inner>
            <inner>two</inner>
        </outer>

    gives you the dict:

    .. code-block:: python

       {u'outer': {u'inner': u'two'}}

    References
    ----------
    * https://github.com/ActiveState/code/tree/master/recipes/Python/116539_turn_structure_XMLdocument
    """
    return node_to_dict(parseString(xmlstring))


def node_to_dict(node):
    """
    Scans through the children of the node and makes a dictionary from the
    content.

    Three cases are differentiated:

    1. If the node contains no other nodes, it is a text-node and
       ``{nodeName: text}`` is merged into the dictionary.
    2. If the node has the attribute ``method`` set to ``true``, then it's children
       will be appended to a list and this list is merged to the dictionary in
       the form: ``{nodeName:list}``.
    3. Else, will call itself recursively on the nodes children
       (merging ``{nodeName: node_to_dict()}`` to the dictionary).

    Parameters
    ----------
    node : `xml.etree.ElementTree.Element`
        A XML element node.

    Returns
    -------
    `dict`:
        The XML element node as a dictionary.
    """
    dic = {}
    for n in node.childNodes:
        if n.nodeType != n.ELEMENT_NODE:
            continue
        if n.getAttribute("multiple") == "true":
            # node with multiple children: put them in a list
            alist = []
            for c in n.childNodes:
                if c.nodeType != n.ELEMENT_NODE:
                    continue
                alist.append(node_to_dict(c))
                dic.update({n.nodeName: alist})
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
    Scans through all children of `~xml.etree.ElementTree.Element` node and
    gathers the text.

    If node has non-text child-nodes, then `~sunpy.util.xml.NotTextNodeError` is raised.

    Parameters
    ----------
    node : `xml.etree.ElementTree.Element`
        A XML element node.

    Returns
    -------
    `str`:
        The `str` context of the XML element node.
    """
    t = ""
    for n in node.childNodes:
        if n.nodeType == n.TEXT_NODE:
            t += n.nodeValue
        else:
            raise NotTextNodeError
    return t
