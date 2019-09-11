# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

import sunpy.util.net


def test_content_disposition_ascii():
    ret = sunpy.util.net.get_content_disposition("Content-Disposition: attachment; filename=foo.txt")
    assert ret == "foo.txt"
    assert isinstance(ret, str)


def test_content_disposition_unicode():
    ret = sunpy.util.net.get_content_disposition("Content-Disposition: attachment; filename*= UTF-8''%e2%82%ac%20rates")
    assert ret == "€ rates"
    assert isinstance(ret, str)


def test_slugify():
    assert sunpy.util.net.slugify("äb c", "b_c")
    assert sunpy.util.net.slugify("file.greg.fits") == "file_greg.fits"
    assert sunpy.util.net.slugify("file.greg.fits", "x") == "filexgreg.fits"
    assert sunpy.util.net.slugify("filegreg.fits") == "filegreg.fits"
    assert sunpy.util.net.slugify("filegreg") == "filegreg"
    assert sunpy.util.net.slugify("f/i*l:e,gr.eg.fits") == "f_i_l_e_gr_eg.fits"
    assert sunpy.util.net.slugify("part1.part2.part3.part4.part5") == "part1_part2_part3_part4.part5"
