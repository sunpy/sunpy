"""
Classes and helper functions for VSO responses.
"""
from collections import defaultdict
from collections.abc import Mapping

import numpy as np
from zeep.helpers import serialize_object

import astropy.units as u

from sunpy.net.base_client import QueryResponseTable
from sunpy.time import parse_time
from sunpy.util._table_attribute import TableAttribute

__all__ = ['VSOQueryResponseTable']


def iter_sort_response(response):
    """
    Sorts the VSO query results by their start time.

    Parameters
    ----------
    response : `zeep.objects.QueryResponse`
        A SOAP Object of a VSO query result

    Returns
    -------
    `list`
        Sorted record items w.r.t. their start time.
    """
    has_time_recs = list()
    has_notime_recs = list()
    for prov_item in response.provideritem:
        if not hasattr(prov_item, 'record') or not prov_item.record:
            continue
        if not hasattr(prov_item.record, 'recorditem') or not prov_item.record.recorditem:
            continue
        rec_item = prov_item.record.recorditem
        for rec in rec_item:
            if hasattr(rec, 'time') and hasattr(rec.time, 'start') and rec.time.start is not None:
                has_time_recs.append(rec)
            else:
                has_notime_recs.append(rec)
    has_time_recs = sorted(has_time_recs, key=lambda x: x.time.start)
    all_recs = has_time_recs + has_notime_recs
    return all_recs


class VSOQueryResponseTable(QueryResponseTable):
    hide_keys = ['fileid', 'fileurl']
    errors = TableAttribute(default=[])

    @classmethod
    def from_zeep_response(cls, response, *, client, _sort=True):
        """
        Construct a table response from the zeep response.
        """
        # _sort is a hack to be able to convert from a legacy QueryResponse to
        # a table response.
        if _sort:
            records = iter_sort_response(response)
        else:
            records = response

        data = []
        for record in records:
            row = defaultdict(lambda: None)
            for key, value in serialize_object(record).items():
                if not isinstance(value, Mapping):
                    if key == "size":
                        # size is in bytes with a very high degree of precision.
                        value = (value * u.Kibyte).to(u.Mibyte).round(5)

                    key = key.capitalize() if key not in cls.hide_keys else key
                    row[key] = value
                else:
                    if key == "wave":
                        # Some records in the VSO have 'kev' which astropy
                        # doesn't recognise as a unit, so fix it.
                        waveunit = value['waveunit']
                        waveunit = 'keV' if waveunit == 'kev' else waveunit

                        row["Wavelength"] = None
                        if value['wavemin'] is not None and value['wavemax'] is not None:
                            row["Wavelength"] = u.Quantity(
                                [float(value['wavemin']), float(value['wavemax'])],
                                unit=waveunit)

                        row["Wavetype"] = value['wavetype']
                        continue
                    for subkey, subvalue in value.items():
                        key_template = f"{key.capitalize()} {subkey.capitalize()}"
                        if key == "time" and subvalue is not None:
                            key_template = f"{subkey.capitalize()} {key.capitalize()}"
                            subvalue = parse_time(subvalue)
                            # Change the display to the 'T'-less version
                            subvalue.format = 'iso'
                        row[key_template] = subvalue
            data.append(row)

        # Reorder the columns to put the most useful ones first.
        data = cls(data, client=client)
        return data._reorder_columns(['Start Time', 'End Time', 'Source',
                                      'Instrument', 'Type', 'Wavelength'],
                                     remove_empty=True)

    def total_size(self):
        if 'size' not in self.colnames:
            return np.nan
        return np.nansum(self['size'])

    def add_error(self, exception):
        self.errors.append(exception)
