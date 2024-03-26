# Author: Florian Mayer <florian.mayer@bitsrc.org>

import sunpy.util.net


def test_content_disposition_ascii():
    ret = sunpy.util.net.get_content_disposition(
        "Content-Disposition: attachment; filename=foo.txt")
    assert ret == "foo.txt"
    assert isinstance(ret, str)


def test_content_disposition_unicode():
    ret = sunpy.util.net.get_content_disposition(
        "Content-Disposition: attachment; filename*= UTF-8''%e2%82%ac%20rates")
    assert ret == "€ rates"
    assert isinstance(ret, str)


def test_slugify():
    assert sunpy.util.net.slugify("ä™") == "äTM"  # Unicode NFKC normalization
    assert sunpy.util.net.slugify("filegreg") == "filegreg"  # no file extension
    assert sunpy.util.net.slugify("filegreg.fits") == "filegreg.fits"  # one file extension
    assert sunpy.util.net.slugify("file.greg.fits") == "file.greg.fits"  # more than one apparent file extension
    assert sunpy.util.net.slugify("AbCdEf") == "AbCdEf"  # uppercase characters
    assert sunpy.util.net.slugify("f/i*l:e,gr.eg.fits") == "f_i_l_e_gr.eg.fits"  # special characters
    assert sunpy.util.net.slugify("file greg'.fits", "x") == "filexgregx.fits"  # custom delimiter
