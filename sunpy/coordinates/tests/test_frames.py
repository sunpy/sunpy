"""
Frame class tests for the coordinates framework.
Modelled after tests in Astropy.
@author: Pritish C. (VaticanCameos)
"""

# NumPy
import numpy as np
from numpy import testing as npt

# Astropy
from astropy import units as u
from astropy.coordinates import (SphericalRepresentation,
                                 CylindricalRepresentation,
                                 CartesianRepresentation)

# Pytest
import pytest

def test_frame_attributes():
    from ..frames import HelioProjective
    from astropy.coordinates.baseframe import FrameAttribute

    class MyProj(HelioProjective):
        # Inherit D0, change d to something else, create a newattr.
        d = (2*u.au).to(u.km)
        newattr = FrameAttribute(default='anattr')

    myproj = MyProj()
    assert myproj.D0 == (1*u.au).to(u.km)
    assert myproj.d == (2*u.au).to(u.km)
    assert myproj.newattr.value == 'anattr'
    assert set(myproj.get_frame_attr_names()) == set([(1*u.au).to(u.km),
                                                      (2*u.au).to(u.km),
                                                      'anattr'])

    myproj = MyProj(D0=1*u.km, d=2*u.km, newattr='changed')
    assert myproj.D0 == 1*u.km
    assert myproj.d == 2*u.km
    assert myproj.newattr.value == 'changed'
    
