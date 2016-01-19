"""
Test Specific Map Scenereos
"""


from __future__ import absolute_import

import os
import pytest
import datetime

import numpy as np

import sunpy.map


def bad_head_test():
    mp=sunpy.map.Map('/Users/mskirk/Desktop/De_Noised_SubPSF_aia_lev1.5_171a_2012_09_23t00_01_59_34z_.fits')
    return mp.units

