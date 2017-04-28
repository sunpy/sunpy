from __future__ import absolute_import

import pytest

from sunpy.net.helio import hec
import sunpy.net.helio.parser as p

def test_suds_unwrapper():
    suds_output = """<?xml version="1.0" encoding="UTF-8"?>
    <S:Envelope ..... >
       <S:Body>
          <helio:queryResponse ... >
             <VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.1" version="1.1">
                <RESOURCE>
                ...
                </RESOURCE>
             </VOTABLE>
          </helio:queryResponse>
       </S:Body>
    </S:Envelope>
    """
    expected_output = """<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.1" version="1.1">
                <RESOURCE>
                ...
                </RESOURCE>
             </VOTABLE>
"""
    assert hec.suds_unwrapper(suds_output) == expected_output

@pytest.mark.online
def test_webservice_parser():
    result = p.webservice_parser()
    assert isinstance(result,list)
