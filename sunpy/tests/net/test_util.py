# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from sunpy.net import util

def test_content_disposition_ascii():
    ret = util.get_content_disposition("Content-Disposition: attachment; filename=foo.txt")
    assert ret == u"foo.txt"
    assert isinstance(ret, unicode)


def test_content_disposition_unicode():
    ret = util.get_content_disposition("Content-Disposition: attachment; filename*= UTF-8''%e2%82%ac%20rates")
    assert ret == u"€ rates"
    assert isinstance(ret, unicode)

def test_slugify():
    assert util.slugify(u"äb c", u"b_c")


def test_slugify_filename():
    assert util.slugify_filename(u"äb c.exe", u"b_c.exe")

def test_slugify_filename2():
    assert util.slugify_filename(u"äb. c.exe", u"b_c.exe")

def test_slugify_filename3():
    assert util.slugify_filename(u"äb. c.exeß", u"b_c.exe")
