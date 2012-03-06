"""
    Provides programs to process and analyze LYRA data. This module is
    currently in development.
   
    Examples
    --------
    To make a LYRA plot
    >>> import sunpy.instr.lyra as lyra
    >>> lobj = lyra.lyra()
    >>> lobj.download('2011/07/23')
    >>> lobj.load()
    >>> lobj.plot()
    
    To Do
    -----
    * An object should be developed to encapsulate all of goes functionality
    * LYRA FITS files store one day's of data at a time.  A simple extension
        of the code would make it easy download multiple days worth of data
        at the one time and create very large arrays of data between
        specific times.
    
    See Also
    --------
    For current data plots see http://proba2.sidc.be
    
    References
    ----------
    | http://proba2.sidc.be/index.html/science/lyra-analysis-manual/article/lyra-analysis-manual?menu=32

    """

import pyfits
import datetime
from matplotlib import pyplot as plt
import urllib,os, copy
from sgmllib import SGMLParser
from sunpy.util.util import anytim
from sunpy.time.util import TimeRange
import numpy as np

"""
class used to scrape the webpage
"""
class URLLister(SGMLParser):
        def reset(self):
                SGMLParser.reset(self)
                self.urls = []

        def start_a(self, attrs):
                href = [v for k, v in attrs if k=='href']
                if href:
                        self.urls.extend(href)

"""
Main LYRA object
"""
class lyra:
    def __init__(self,time):
        self.verbose = False
        self.filename = None
        self.data = None
        self.columns = None
        self.level = 2
        self.dataType = 'std'
        self.downloadto = os.path.expanduser('~')+os.sep
        self.location = 'http://proba2.oma.be/lyra/data/bsd/'
        self.prefix = 'lyra_'
        self.time = anytim(time)
        self.tstart = self.time
        self.tend = self.tstart + datetime.timedelta(days = 1)
        self.nt = None
        
    def download(self):
        """ Function to download LYRA data, and/or set the filename where it can be found"""

        self.time = anytim(self.time)

        # date-based subdirectory
        dateLocation = self.time.strftime('%Y/%m/%d/')

        # date-based filename
        dateFilename = self.time.strftime('%Y%m%d-')

        # extension to the file name
        extension = '000000_lev'+ str(self.level) + '_' + self.dataType + '.fits'

        # calculated file name
        requestedLocation = self.location + dateLocation
        requestedFile  = self.prefix + dateFilename + extension

        f = os.path.expanduser(self.downloadto) + os.sep + requestedFile

        isAlreadyThere = os.path.isfile(f)

        if isAlreadyThere:
            if self.verbose:
                print('File '+ f + ' already exists.')
            self.filename = f
        else:
            self.filename = wgetDownload(requestedLocation,requestedFile, self.downloadto)

    def load(self):
        """Load the FITS file data into the object"""
        hdu = pyfits.open(self.filename)
        self.data = hdu[1].data
        self.columns = hdu[1].columns
        hdu.close()

        # define the start and end times 
        self.tstart = self.time
        self.tend = self.time + datetime.timedelta(seconds = self.data.field(0)[-1])
        self.nt = self.data.field(0).size

    def plot(self):
        """Plot the LYRA data"""
        names = self.columns.names
        units = self.columns.units
        time = self.data.field(0)

        fig, axs = plt.subplots(4,1, sharex=True)

        for j in range(1,5):
            axs[j-1].plot(time,self.data.field(j))
            axs[j-1].set_ylabel( names[j] + ' (' + units[j] + ')' )

        axs[0].set_title(str(self.tstart) + ' - '+ str(self.tend))
        plt.xlabel( names[0] + ' (' + units[0] + ')' )
        plt.show()

#
# general purpose downloader using wget
#
def wgetDownload(requestedLocation,requestedFile, downloadto, verbose = False, overwrite=False):

    usock = urllib.urlopen(requestedLocation)
    parser = URLLister()
    parser.feed(usock.read())
    usock.close()
    parser.close()

    if verbose:
        print('Requested file = '+requestedFile)

    print(parser.urls)

    if not requestedFile in parser.urls:
        if verbose:
            print('Requested file not found.')
        return None
    else:
        if verbose:
            print('Downloading '+ requestedLocation+requestedFile)
        remoteBaseURL = '-B ' + requestedLocation + ' '
        localDir = ' -P'+downloadto + ' '
        command = 'wget -r -l1 -nd --no-parent ' + requestedLocation+requestedFile + localDir + remoteBaseURL
        if verbose:
            print('Executing '+command)
            print('File located at ' + downloadto + requestedFile)
        os.system(command)
        return downloadto + requestedFile


def sub(lobj,t1,t2):
    """Returns a subsection of LYRA data based on times
         which is also a LYRA object."""
    tr = TimeRange(t1,t2)
    if tr.t1 < lobj.tstart:
        print('Requested start time less than data start time.')
        indexStart = 0
    else:
        indexStart = getLyraIndexGivenTime(lobj, tr.t1)
    if tr.t2 > lobj.tend:
        print('Requested end time greater than data start time.')
        indexEnd = lobj.nt-1
    else:
        indexEnd = getLyraIndexGivenTime(lobj, tr.t2)
    mask1 = lobj.data.field(0) > lobj.data.field(0)[indexStart]
    mask2 = lobj.data.field(0) < lobj.data.field(0)[indexEnd]
    newlobj = lyra(lobj.time)
    newlobj.tstart = getLyraTimeGivenIndex(lobj,indexStart)
    newlobj.tend = getLyraTimeGivenIndex(lobj,indexEnd)
    newlobj.data = lobj.data[mask1*mask2]
    newlobj.nt = newlobj.data.field(0).size
    firstTime = newlobj.data.field(0)[0]
    newlobj.data.field(0)[:] = newlobj.data.field(0)[:] - firstTime
    newlobj.columns = lobj.columns
    return newlobj


def getLyraTimeGivenIndex(lobj,index):
    """Given an index, return a time from the LYRA sample times"""
    return lobj.time + datetime.timedelta(seconds = lobj.data.field(0)[index])

def getLyraIndexGivenTime(lobj,targetTime):
    """Given a target time, return an index based on the LYRA sample times"""
    earlyIndex = 0
    lateIndex = lobj.nt-1
    midIndex = 0.5*(earlyIndex+lateIndex)
    while ( (lateIndex-earlyIndex) > 1):
        midIndex = np.around(0.5*(earlyIndex + lateIndex))
        tmidpoint = getLyraTimeGivenIndex(lobj,midIndex)
        if targetTime >= tmidpoint:
            earlyIndex = midIndex
        if targetTime < tmidpoint:
            lateIndex = midIndex
    return midIndex


def sumup(lobj,binsize, average = True):
    """Add up the data using a bin size in seconds"""
    mask = np.zeros(shape=(lobj.nt,),dtype='bool')
    binStart = np.array([],dtype='int')
    binEnd = np.array([],dtype='int')
    mask[:] = False
    mask[0] = True
    copydata = copy.deepcopy(lobj.data)
    newlobj = lyra(lobj.time)
    newlobj.columns = lobj.columns
    iStart = 0
    nNew = 0
    # Count the number of elements in the summed array
    while iStart < lobj.nt-1:    
        iCount = sumupFind(lobj,binsize,iStart)
        binStart = np.append(binStart,[iStart])
        binEnd = np.append(binEnd,[iStart+iCount])
        print binStart, binEnd
        mask[iStart + iCount] = True
        copydata.field(1)[iStart] = np.average(lobj.data.field(1)[ binStart[nNew]:binEnd[nNew] ] )
        copydata.field(2)[iStart] = np.average(lobj.data.field(2)[ binStart[nNew]:binEnd[nNew] ] )
        copydata.field(3)[iStart] = np.average(lobj.data.field(3)[ binStart[nNew]:binEnd[nNew] ] )
        copydata.field(4)[iStart] = np.average(lobj.data.field(4)[ binStart[nNew]:binEnd[nNew] ] )
        iStart = iStart + iCount + 1
        nNew = nNew + 1
    print mask

    newlobj = lyra(lobj.time)
    newlobj.tstart = getLyraTimeGivenIndex(lobj,0)
    newlobj.tend = getLyraTimeGivenIndex(lobj,1)
    newlobj.data = copydata[binStart[:]]
    newlobj.nt = newlobj.data.field(0).size
    firstTime = newlobj.data.field(0)[0]
    newlobj.data.field(0)[:] = newlobj.data.field(0)[:] - firstTime
    newlobj.columns = lobj.columns

    return newlobj

def sumupFind(lobj,binsize,iStart):
    i = 0
    iStartTime = getLyraTimeGivenIndex(lobj,iStart)
    while (iStart + i < lobj.nt-1) and (getLyraTimeGivenIndex(lobj,iStart+i) < iStartTime + datetime.timedelta(seconds = binsize)):
        i = i + 1
    return i-1

def subthensum(lobj,t1,t2,binsize,average = True):
    """Subsection the time series then sum it up"""
    n1 = sub(lobj,t1,t2)
    return sumup(n1,binsize,average=average)
    