from datetime import datetime

from dateutil.relativedelta import relativedelta

from sunpy.extern.parse import parse
from sunpy.net.scraper import get_timerange_from_exdict


def check_timerange(pattern, url, timerange):
    """
    Checks whether the time range represented in *url* intersects
    with the given time range.

    Parameters
    ----------
    url : `str`
        URL of the file.
    timerange : `~sunpy.time.TimeRange`
        Time interval for which files were searched.

    Returns
    -------
    `bool`
        `True` if URL's valid time range overlaps the given timerange, else `False`.
    """
    exdict = parse(pattern, url).named
    if exdict['year'] < 100:
        exdict['year'] = 2000 + exdict['year']
    if 'month' not in exdict:
                if 'month_name' in exdict:
                    exdict['month'] = datetime.strptime(exdict['month_name'], '%B').month
                elif 'month_name_abbr' in exdict:
                    exdict['month'] = datetime.strptime(exdict['month_name_abbr'], '%b').month
    tr = get_timerange_from_exdict(exdict)
    return tr.intersects(timerange)

def smallerPattern(directoryPattern):
    """
    Obtain the smaller time step for the given pattern.
    """
    try:
        if "%S" in directoryPattern:
            return relativedelta(seconds=1)
        elif "%M" in directoryPattern:
            return relativedelta(minutes=1)
        elif any(hour in directoryPattern for hour in ["%H"]):
            return relativedelta(hours=1)
        elif any(day in directoryPattern for day in ["%d", "%j"]):
            return relativedelta(days=1)
        elif any(month in directoryPattern for month in ["%b", "%B", "%m"]):
            return relativedelta(months=1)
        elif any(year in directoryPattern for year in ["%Y", "%y"]):
            return relativedelta(years=1)
        else:
            return None
    except Exception:
        raise
