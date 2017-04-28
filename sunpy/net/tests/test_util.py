# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

import sunpy.util.net
from sunpy.extern import six

def test_content_disposition_ascii():
    ret = sunpy.util.net.get_content_disposition("Content-Disposition: attachment; filename=foo.txt")
    assert ret == u"foo.txt"
    assert isinstance(ret, six.text_type)


def test_content_disposition_unicode():
    ret = sunpy.util.net.get_content_disposition("Content-Disposition: attachment; filename*= UTF-8''%e2%82%ac%20rates")
    assert ret == u"€ rates"
    assert isinstance(ret, six.text_type)

def test_slugify():
    assert sunpy.util.net.slugify(u"äb c", u"b_c")
    assert sunpy.util.net.slugify(u"file.greg.fits") == u"file_greg.fits"
    assert sunpy.util.net.slugify(u"file.greg.fits", u"x") == u"filexgreg.fits"
    assert sunpy.util.net.slugify(u"filegreg.fits") == u"filegreg.fits"
    assert sunpy.util.net.slugify(u"filegreg") == u"filegreg"
    assert sunpy.util.net.slugify(u"f/i*l:e,gr.eg.fits") == u"f_i_l_e_gr_eg.fits"
    assert sunpy.util.net.slugify(u"part1.part2.part3.part4.part5") == u"part1_part2_part3_part4.part5"
