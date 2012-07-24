# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

import os
import shutil

from tempfile import mkdtemp
from datetime import datetime

import numpy as np

from sunpy.data.sample import CALLISTO_IMAGE
from sunpy.spectra.sources.callisto import (
    CallistoSpectrogram, query, download
)


def test_read():
    ca = CallistoSpectrogram.read(CALLISTO_IMAGE)
    assert ca.start == datetime(2011, 9, 22, 10, 30, 0, 51000)
    assert (
        ca.t_init ==
        (datetime(2011, 9, 22, 10, 30) - datetime(2011, 9, 22)).seconds
    )
    assert ca.shape == (200, 3600)
    assert ca.t_delt == 0.25
    # Test linearity of time axis.
    assert np.array_equal(
        ca.time_axis, np.linspace(0, 0.25 * (ca.shape[1] - 1), ca.shape[1])
    )
    assert ca.dtype == np.uint8


def test_query():
    URL = 'http://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto/2011/09/22/'

    result = list(query(
        datetime(2011, 9, 22, 5), datetime(2011, 9, 22, 6), set(["BIR"])
    ))
    RESULTS = [
        "BIR_20110922_050000_01.fit.gz",
        "BIR_20110922_051500_01.fit.gz",
        "BIR_20110922_053000_01.fit.gz",
        "BIR_20110922_060000_01.fit.gz",
        "BIR_20110922_050000_03.fit.gz",
        "BIR_20110922_051500_03.fit.gz",
        "BIR_20110922_053000_03.fit.gz",
        "BIR_20110922_054500_03.fit.gz",
        "BIR_20110922_060000_03.fit.gz"
    ]

    RESULTS.sort()

    assert result == [URL + res for res in RESULTS]


def test_query_number():
    URL = 'http://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto/2011/09/22/'

    result = list(query(
        datetime(2011, 9, 22, 5), datetime(2011, 9, 22, 6), set(["BIR"]), 1
    ))
    RESULTS = [
        "BIR_20110922_050000_01.fit.gz",
        "BIR_20110922_051500_01.fit.gz",
        "BIR_20110922_053000_01.fit.gz",
        "BIR_20110922_060000_01.fit.gz",
    ]

    RESULTS.sort()

    assert result == [URL + res for res in RESULTS]


def test_download():
    directory = mkdtemp()
    try:
        result = query(
            datetime(2011, 9, 22, 5), datetime(2011, 9, 22, 6), set(["BIR"]), 1
        )
        RESULTS = [
            "BIR_20110922_050000_01.fit.gz",
            "BIR_20110922_051500_01.fit.gz",
            "BIR_20110922_053000_01.fit.gz",
            "BIR_20110922_060000_01.fit.gz",
        ]
        download(result, directory)
        assert sorted(os.listdir(directory)) == RESULTS
    finally:
        shutil.rmtree(directory)
