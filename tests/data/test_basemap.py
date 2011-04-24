"""BaseMap tests"""
import unittest
import sunpy
import pyfits

class TestBaseMap(unittest.TestCase):
    """Tests the BaseMap class"""
    def setUp(self):
        self.file = 'doc/sample-data/AIA20110319_105400_0171.fits'
        self.map = sunpy.Map(self.file) 

    def tearDown(self):
        self.map = None

    def test_fits_data_comparison(self):
        """Make sure the data is the same in pyfits and SunPy"""
        fits = pyfits.open(self.file)
        self.assertEqual(self.map.tolist(), fits[0].data.tolist(),
                         'data not preserved')
    def test_fits_header_comparison(self):
        """Make sure the header is the same in pyfits and SunPy"""
        fits = pyfits.open(self.file)
        
        # Access fits data once to apply scaling-related changes and update
        # header information in fits[0].header
        fits[0].data

        self.assertEqual(self.map.header.keys(), fits[0].header.keys(),
                         'header not preserved')

suite = unittest.TestLoader().loadTestsFromTestCase(TestBaseMap)
unittest.TextTestRunner(verbosity=2).run(suite)
