from datetime import timedelta

import astropy.units as u

from sunpy.net.base_client import QueryResponseTable
from sunpy.time import parse_time


def create_norh_url(date, wavelength):
    year = date.strftime("%Y")
    month = date.strftime("%m")
    day = date.strftime("%d")
    if wavelength == 34 * u.GHz:
        freq = 'tcz'
    elif wavelength == 17 * u.GHz:
        freq = 'tca'
    value_ = year[-2:] + month + day
    base = f"https://solar.nro.nao.ac.jp/norh/data/tcx/{year}/{month}/{freq}{value_}"
    return (base)


def mock_query_object(timerange, client):
    """
    Creating a Query Response object and prefilling it with some information
    """
    # Creating a Query Response Object
    start = timerange.start
    end = timerange.end
    wave = 17*u.GHz
    delta = end-start
    resp = []
    for i in range(delta.datetime.days + 1):
        start_time = start.datetime + timedelta(days=i)
        end_time = start_time + timedelta(days=1) - timedelta(milliseconds=1)
        obj = {
            'Start Time': parse_time(start_time),
            'End Time': parse_time(end_time),
            'Instrument': 'NORH',
            'Source': 'NAOJ',
            'Provider': 'NRO',
            'Wavelength': wave,
            'url': create_norh_url(start_time, wave)
        }
        resp.append(obj)
    results = QueryResponseTable(resp)
    results.client = client
    return results
