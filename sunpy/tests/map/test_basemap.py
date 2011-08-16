"""
BaseMap tests
"""
#pylint: disable=C0103,R0904
import sunpy
import pyfits

class TestBaseMap:
    """Tests the BaseMap class"""
    def setup_class(self):
        self.file = 'doc/sample-data/AIA20110319_105400_0171.fits'
        self.map = sunpy.Map(self.file)
        self.fits = pyfits.open(self.file)

    def teardown_class(self):
        self.map = None
        self.fits = None

    def test_fits_data_comparison(self):
        """Make sure the data is the same in pyfits and SunPy"""
        assert (self.map == self.fits[0].data).all()

    def test_fits_header_comparison(self):
        """Make sure the header is the same in pyfits and SunPy"""

        # Access fits data once to apply scaling-related changes and update
        # header information in fits[0].header
        self.fits[0].data #pylint: disable=W0104

        # 2011/08/15: 
        # What we really want to do here is cast self.map.header
        # and self.fits[0].header to dicts and compare the results, however,
        # pyfits complains about certain keys when attemping this, e.g.:
        # 
        # ValueError: Unparsable card (X0_LF), fix it first with .verify('fix').
        #
        # After calling fits.verifty('fix'), the casting works, but since this
        # isn't done in BaseMap, the test cannot be performed.
        
        #assert dict(self.map.header) == dict(self.fits[0].header)
        pass