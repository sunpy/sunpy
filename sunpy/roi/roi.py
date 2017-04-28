from __future__ import absolute_import

from sunpy.time import TimeRange
from sunpy.time import parse_time

__all__ = ['roi']


class roi(object):
    """
    A generalized region of interest (ROI) object

    Parameters
    ----------
    times : list (optional)
        A list of 1 or 2 parse_time-readable times
    description : str (optional)
        A text description of the ROI
    source : str (optional)
        A description of where this ROI comes from
        (e.g. the instrument, 'RHESSI', 'LYRA LYTAF')

    Attributes
    ----------
    start_time : datetime object containing the start time of the ROI
    end_time : datetime object containing the end time of the ROI
    description : A string descriptor of the ROI event type
        (e.g. 'attenuator change', 'LAR', 'SAA', 'flare')
    source : A string descriptor of the ROI source (e.g. 'LYRA', 'RHESSI')


    Methods
    -------
    time_range()
        Return a time range object from the start and end times of the ROI

    Examples
    --------
    >>> from sunpy.roi import roi

    >>> result = roi(times=['2011-02-15 04:34:09','2011-02-15 04:48:21'], description='UV occult.',source='LYRA LYTAF')

    >>> result = roi(times='2013-05-12 03:12:00')
    """

    def __init__(self, times=None, description=None, source=None):
        # time could be a list with one or two elements
        if times and type(times) == list:
            if len(times) == 1:
                # if only one time given, make start and end times the same
                self.start_time = parse_time(times[0])
                self.end_time = parse_time(times[0])
            elif len(times) == 2:
                self.start_time = parse_time(times[0])
                self.end_time = parse_time(times[1])
            else:
                self.start_time = None
                self.end_time = None
        elif type(times) == str:
            self.start_time = parse_time(times)
            self.end_time = parse_time(times)
        else:
            self.start_time = None
            self.end_time = None

        # description of the ROI event type
        if description:
            self.description = str(description)
        else:
            self.description = None

        # optional description of where the ROI came from
        if source == None:
            self.source = "Unknown"
        else:
            self.source = source

    def time_range(self):
        """Returns a TimeRange using the start and end times"""
        if self.start_time and self.end_time:
            return TimeRange(self.start_time, self.end_time)

    def __repr__(self):
        """Print out info on the ROI"""
        if not self.start_time:
            startstring = 'None'
        else:
            startstring = self.start_time.isoformat()

        if not self.end_time:
            endstring = 'None'
        else:
            endstring = self.end_time.isoformat()
        return('SunPy Region-of-interest (ROI) object' +
        '\n-------------------------------------' +
        '\nSource: \t\t' + self.source +
        '\nStart time:\t\t' + startstring +
        '\nEnd time: \t\t' + endstring +
        '\nEvent description:\t' + str(self.description))
