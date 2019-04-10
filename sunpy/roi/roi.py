"""
This module implements a basic region of interest class.
"""
import textwrap

from sunpy.time import TimeRange, parse_time

__all__ = ['roi']


class roi(object):
    """
    A generalized Region Of Interest (ROI) class.

    Parameters
    ----------
    times : `list`, optional
        A list of 1 or 2 `sunpy.time.parse_time` readable times.
    description : `str`, optional
        A text description of the ROI.
    source : `str`, optional
        A description of of the source where this ROI comes from.

    Attributes
    ----------
    start_time : `astropy.time.Time`
        The start time of the ROI.
    end_time : `astropy.time.Time`
        The end time of the ROI.
    description : `str`
        A string descriptor of the ROI event type (e.g., 'attenuator change', 'LAR', 'SAA', 'flare').
    source : `str`
        A string descriptor of the ROI source (e.g., 'LYRA', 'RHESSI').

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
        if source is None:
            self.source = "Unknown"
        else:
            self.source = source

    def time_range(self):
        """
        Returns a `sunpy.time.TimeRange` using the start and end times.
        """
        if self.start_time and self.end_time:
            return TimeRange(self.start_time, self.end_time)

    def __repr__(self):
        """
        Print out info on the ROI.
        """
        if not self.start_time:
            startstring = 'None'
        else:
            startstring = self.start_time.iso

        if not self.end_time:
            endstring = 'None'
        else:
            endstring = self.end_time.iso
        return textwrap.dedent(f"""\
        SunPy Region-of-interest (ROI) object
        -------------------------------------
        Source:            {self.source}
        Start time:        {startstring}
        End time:          {endstring}
        Event description: {self.description}
        """)
