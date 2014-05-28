from sunpy.util import unit_conversion as ut
from numpy.testing import assert_array_almost_equal as a
from astropy import units as u

def test_kelvin_to_keV():
    a(ut.kelvin_to_keV(300 * u.K), 2.585e-05 * u.keV, decimal = 7)
    a(ut.kelvin_to_keV(273.5 * u.K), 2.356e-05 * u.keV, decimal = 7)
    a(ut.kelvin_to_keV(3023.12 * u.K), 2.605e-04 * u.keV, decimal = 7) 


	
	
