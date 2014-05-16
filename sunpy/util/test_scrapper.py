import pytest
import datetime
import os

import sunpy.data.test
from sunpy.util.scrapper import Scrapper

PATTERN_EXAMPLES = [
    ('%b%y', datetime.timedelta(days=31)),
    ('%m%y', datetime.timedelta(days=31)),
    ('%H%d', datetime.timedelta(hours=1)),
    ]

class TestScrapper:
    def testDirectoryDatePattern(self):
        s = Scrapper('%Y/%m/%d/%Y%m%d_%H%M%S_59.fit.gz')
        testpath = '2014/03/05/20140305_013000_59.fit.gz'
        d = datetime.datetime(2014,3,5,1,30)
        assert s.matches(testpath, d)
        
    def testDirectoryDatePatternFalse(self):
        s = Scrapper('%Y/%m/%d/%Y%m%d_%H%M%S_59.fit.gz')
        testpath = '2013/03/05/20140305_013000_59.fit.gz'
        d = datetime.datetime(2014,3,5,1,30)
        assert not s.matches(testpath, d)
        
    def testDirectoryObsPattern(self):
        s = Scrapper('%y%m%d/{observatory}_%Y%m%d.fits', observatory = 'SDO')
        testpath = '140305/SDO_20140305.fits'
        d = datetime.datetime(2014,3,5)
        assert s.matches(testpath, d)
        
    def testDirectoryRange(self):
        s = Scrapper('%Y/%m/%d/%Y%m%d_%H.fit.gz')
        directory_list = ['2009/12/30/', '2009/12/31/', '2010/01/01/', 
                          '2010/01/02/', '2010/01/03/']
        startdate = datetime.datetime(2009,12,30)
        enddate = datetime.datetime(2010,1,3)
        assert s.range(startdate, enddate) == directory_list
        
    def testDirectoryRangeFalse(self):
        s = Scrapper('%Y%m%d/%Y%m%d_%H.fit.gz')
        directory_list = ['20091230/', '20091231/', '20100101/', 
                          '20090102/', '20090103/']
        startdate = datetime.datetime(2009,12,30)
        enddate = datetime.datetime(2010,1,3)
        assert s.range(startdate, enddate) != directory_list

    @pytest.mark.parametrize('pattern, mintime', PATTERN_EXAMPLES)
    def test_smallerPattern(self, pattern, mintime):
        assert  mintime == Scrapper('')._smallerPattern(pattern)

    def testDirectoryRangeHours(self):
        s = Scrapper('%Y%m%d_%H/%H%M.csv')
        startdate = datetime.datetime(2009,12,31,23,40)
        enddate = datetime.datetime(2010,01,01,01,15)
        assert len(s.range(startdate, enddate)) == 3 #3 directories (1 per hour)

    def testDirectoryRange_single(self):
        s = Scrapper('%Y%m%d/%H_%M.csv')
        startdate = datetime.datetime(2010,10,10,5,00)
        enddate = datetime.datetime(2010,10,10,7,00)
        assert len(s.range(startdate, enddate)) == 1
        
    def testDirectoryRange_Month(self):
        s = Scrapper('%Y%m/%d/%j_%H.txt')
        startdate = datetime.datetime(2008, 2,20,10)
        enddate = datetime.datetime(2008, 3, 2, 5)
        assert len(s.range(startdate, enddate)) == 12
        startdate = datetime.datetime(2009, 2,20,10)
        enddate = datetime.datetime(2009, 3, 2, 5)
        assert len(s.range(startdate, enddate)) == 11

    def testNoDirectory(self):
        s = Scrapper('files/%Y%m%d_%H%M.dat')
        startdate = datetime.datetime(2010, 1, 10, 20, 30)
        enddate = datetime.datetime(2010, 1, 20, 20, 30)
        assert len(s.range(startdate, enddate)) == 1

    def testExtractDates_usingPattern(self):
        s = Scrapper('data/%Y/%m/%d/fits/swap/swap_00174_fd_%Y%m%d_%H%M%S.fts.gz')
        testURL = 'data/2014/05/14/fits/swap/swap_00174_fd_20140514_200135.fts.gz'
        timeURL = datetime.datetime(2014, 5, 14, 20, 1, 35)
        assert s.extractDateURL(testURL) == timeURL

    def testURL_pattern(self):
        s = Scrapper('fd_%Y%m%d_%H%M%S.fts')
        assert s.URL_followsPattern('fd_20130410_231211.fts') == True
        assert s.URL_followsPattern('fd_20130410_231211.fts.gz') == False
        assert s.URL_followsPattern('fd_20130410_ar_231211.fts.gz') == False
        
    # def testFilesRange_sameDirectory_local(self):
    #     s = Scrapper(os.path.join('file://',sunpy.data.test.rootdir,
    #                               'EIT','efz%Y%m%d.%H%M%S_s.fits'))
    #     startdate = datetime.datetime(2004, 3, 1, 4, 0)
    #     enddate = datetime.datetime(2004, 3, 1, 6, 30)
    #     assert len(s.filelist(startdate, enddate)) == 3
    #     startdate = datetime.datetime(2010, 1, 10, 20, 30)
    #     enddate = datetime.datetime(2010, 1, 20, 20, 30)
    #     assert len(s.filelist(startdate, enddate)) == 0

    def testFilesRange_sameDirectory_remote(self):
        s = Scrapper('http://solarmonitor.org/data/%Y/%m/%d/fits/{instrument}/{instrument}_00174_fd_%Y%m%d_%H%M%S.fts.gz', instrument = 'swap')
        startdate = datetime.datetime(2014, 5, 14, 0, 0)
        enddate = datetime.datetime(2014, 5, 14, 6, 30)
        assert len(s.filelist(startdate, enddate)) == 2
        startdate = datetime.datetime(2014, 5, 14, 21, 0)
        enddate = datetime.datetime(2014, 5, 14, 23, 30)
        assert len(s.filelist(startdate, enddate)) == 0


    

               


        #TODO: data gaps
                        

#def main():
#    pytest.main()

#if __name__ == '__main__':
#    main()
