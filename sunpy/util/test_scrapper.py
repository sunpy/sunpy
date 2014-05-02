import pytest
import datetime
from .scrapper import Scrapper

pattern_examples = [
    ('%b%y', datetime.timedelta(days=31)),
    ('%m%y', datetime.timedelta(days=31)),
    ]

class ScrapperTests:

    def testDirectoryDatePattern(self):
        s = Scrapper('%Y/%m/%d/%Y%m%d_%H%M%S_59.fit.gz')
        testpath = '2014/03/05/20140305_013000_59.fit.gz'
        d = datetime.datetime(2014,3,5,1,30)
        self.failUnless(s.matches(testpath, d))
    def testDirectoryDatePatternFalse(self):
        s = Scrapper('%Y/%m/%d/%Y%m%d_%H%M%S_59.fit.gz')
        testpath = '2013/03/05/20140305_013000_59.fit.gz'
        d = datetime.datetime(2014,3,5,1,30)
        self.failIf(s.matches(testpath, d))
    def testDirectoryObsPattern(self):
        s = Scrapper('%y%m%d/{observatory}_%Y%m%d.fits', observatory = 'SDO')
        testpath = '140305/SDO_20140305.fits'
        d = datetime.datetime(2014,3,5)
        self.failUnless(s.matches(testpath, d))
    def testDirectoryRange(self):
        s = Scrapper('%Y/%m/%d/%Y%m%d_%H.fit.gz')
        directory_list = ['2009/12/30/', '2009/12/31/', '2010/01/01/', 
                          '2010/01/02/', '2010/01/03/']
        startdate = datetime.datetime(2009,12,30)
        enddate = datetime.datetime(2010,1,3)
        self.failUnless(s.range(startdate, enddate) == directory_list)
    def testDirectoryRangeFalse(self):
        s = Scrapper('%Y%m%d/%Y%m%d_%H.fit.gz')
        directory_list = ['20091230/', '20091231/', '20100101/', 
                          '20090102/', '20090103/']
        startdate = datetime.datetime(2009,12,30)
        enddate = datetime.datetime(2010,1,3)
        self.failIf(s.range(startdate, enddate) == directory_list)
    @pytest.mark.parametrize('pattern, mintime', pattern_examples)
    def test_smallerPattern(self):
        assert  mintime == Scrapper('')._smallerPattern(pattern)
     
# check when there's not directory structure (no /s)         
                        

def main():
    unittest.main()

if __name__ == '__main__':
    main()
