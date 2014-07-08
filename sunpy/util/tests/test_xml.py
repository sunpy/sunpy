from __future__ import absolute_import

from xml.parsers.expat import ExpatError
from xml.dom.minidom import Document
from xml.dom.minidom import parseString

import pytest

from sunpy.util import xml


def test_xml_to_dict1():
    """
    should return dict of xml string
    """
    source_xml = "<outer>\
          <inner1>one</inner1>\
          <inner2>two</inner2>\
        </outer>"

    xml_dict = xml.xml_to_dict(source_xml)
    expected_dict = {u'outer': {u'inner2': u'two', u'inner1': u'one'}}

    assert xml_dict == expected_dict

def test_xml_to_dict2():
    """
    should return dict of xml string
    and if a tag is duplicated it takes the last one.
    """
    source_xml = "<outer>\
                    <inner1>one-one</inner1>\
                    <inner1>one-two</inner1>\
                    <inner2>two-one</inner2>\
                    <inner2>two-two</inner2>\
                 </outer>"

    xml_dict = xml.xml_to_dict(source_xml)
    expected_dict = {u'outer': {u'inner2': u'two-two', u'inner1': u'one-two'}}

    assert xml_dict == expected_dict

def test_xml_to_dict3():
    """
    should return dict of xml string
    with empty value if there are no inner elements
    """
    source_xml = "<outer/>"

    xml_dict = xml.xml_to_dict(source_xml)
    expected_dict = {u'outer': ''}

    assert xml_dict == expected_dict

def test_xml_to_dict4():
    """
    should return dict of xml string
    with empty value if there are no inner elements
    """
    source_xml = "<outer></outer>"

    xml_dict = xml.xml_to_dict(source_xml)
    expected_dict = {u'outer': ''}

    assert xml_dict == expected_dict

def test_xml_to_dict5():
    """
    should return dict of xml string
    with 2 layer nesting
    """
    source_xml = "<outer>\
                    <mid1>\
                        <inner1>one-one</inner1>\
                    </mid1>\
                    <mid2>\
                        <inner2>two-one</inner2>\
                    </mid2>\
                 </outer>"

    xml_dict = xml.xml_to_dict(source_xml)
    expected_dict = {u'outer': {u'mid2': {u'inner2': u'two-one'}, u'mid1': {u'inner1': u'one-one'}}}

    assert xml_dict == expected_dict

def test_xml_to_dict6():
    """
    should return dict of xml string
    with 2 layer nesting and if a tag is duplicated it takes the last one.
    """
    source_xml = "<outer>\
                    <mid>\
                        <inner1>one-one</inner1>\
                    </mid>\
                    <mid>\
                        <inner2>two-one</inner2>\
                    </mid>\
                 </outer>"

    xml_dict = xml.xml_to_dict(source_xml)
    expected_dict = {u'outer': {u'mid': {u'inner2': u'two-one'}}}

    assert xml_dict == expected_dict

def test_xml_to_dict7():
    """
    should raise TypeError when passed None
    """
    assert pytest.raises(TypeError, xml.xml_to_dict, None)

def test_xml_to_dict8():
    """
    should raise TypeError when passed non string
    """
    assert pytest.raises(TypeError, xml.xml_to_dict, 9)

def test_xml_to_dict9():
    """
    should raise ExpatError when passed empty string
    """
    assert pytest.raises(ExpatError, xml.xml_to_dict, "")

def test_xml_to_dict10():
    """
    should raise ExpatError when passed space
    """
    assert pytest.raises(ExpatError, xml.xml_to_dict, " ")


def test_get_node_text1():
    """
    should raise NotTextNodeError if there is a non text node.
    """
    doc = Document()
    outer = doc.createElement("outer")
    doc.appendChild(outer)
    pytest.raises(xml.NotTextNodeError, xml.get_node_text, doc)

def test_get_node_text2():
    """
    should return empty string for a node with no child nodes.
    """
    assert xml.get_node_text(Document()) == ""

def test_get_node_text3():
    """
    should return node text
    """
    node = parseString("<outer>one</outer>")
    text_node = node.childNodes[0]

    assert xml.get_node_text(text_node) == "one"

def test_get_node_text4():
    """
     should raise AttributeError when sent None
    """
    assert pytest.raises(AttributeError, xml.get_node_text, None)

def test_get_node_text5():
    """
     should raise AttributeError when sent wrong type
    """
    assert pytest.raises(AttributeError, xml.get_node_text, "wrong type")

def test_node_to_dict1():
    """
    should return dict of node
    """

    doc = Document()

    outer = doc.createElement("outer")
    doc.appendChild(outer)

    inner1 = doc.createElement("inner1")
    inner2 = doc.createElement("inner2")
    outer.appendChild(inner1)
    outer.appendChild(inner2)

    inner1_text = doc.createTextNode("one")
    inner2_text = doc.createTextNode("two")
    inner1.appendChild(inner1_text)
    inner2.appendChild(inner2_text)

    expected_dict = {'outer': {'inner2': 'two', 'inner1': 'one'}}
    xml_dict = xml.node_to_dict(doc)

    assert xml_dict == expected_dict

def test_node_to_dict2():
    """
    should return dict of node double nested
    """

    doc = Document()

    outer = doc.createElement("outer")
    doc.appendChild(outer)

    mid1 = doc.createElement("mid1")
    outer.appendChild(mid1)
    mid2 = doc.createElement("mid2")
    outer.appendChild(mid2)

    inner1 = doc.createElement("inner1")
    inner2 = doc.createElement("inner2")
    mid1.appendChild(inner1)
    mid2.appendChild(inner2)

    inner1_text = doc.createTextNode("one")
    inner2_text = doc.createTextNode("two")
    inner1.appendChild(inner1_text)
    inner2.appendChild(inner2_text)

    expected_dict = {'outer': {'mid2': {'inner2': 'two'}, 'mid1': {'inner1': 'one'}}}
    xml_dict = xml.node_to_dict(doc)

    assert xml_dict == expected_dict

def test_node_to_dict3():
    """
    should return empty dict when sent empty doc
    """
    expected_dict = {}
    xml_dict = xml.node_to_dict(Document())

    assert xml_dict == expected_dict

def test_node_to_dict4():
    """
    should raise AttributeError when sent wrong type
    """
    assert pytest.raises(AttributeError, xml.node_to_dict, 9)

def test_node_to_dict5():
    """
    should raise AttributeError when sent None
    """
    assert pytest.raises(AttributeError, xml.node_to_dict, None)
