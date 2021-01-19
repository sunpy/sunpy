"""
Classes and helper functions for VSO responses.
"""
from collections import defaultdict

from zeep.helpers import serialize_object

import astropy.units as u

from sunpy.net.base_client import QueryResponseTable
from sunpy.time import parse_time
from sunpy.util.table_attribute import TableAttribute

__all__ = ['VSOQueryResponseTable']


def iter_sort_response(response):
    """
    Sorts the VSO queryresults by their start time.

    Parameters
    ----------
    response : `zeep.objects.QueryResponse`
        A SOAP Object of a VSO queryresult

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

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # If we haven't been given a client create one
        if self.client is None:
            # Import here to avoid circular import
            from sunpy.net.vso import VSOClient
            self.client = VSOClient()

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
                if not isinstance(value, dict):
                    if key == "size":
                        # size is in bytes with a very high degree of precision.
                        value = (value * u.Kibyte).to(u.Mibyte).round(5)

                    key = key.capitalize() if key in cls.hide_keys else key
                    row[key] = value
                else:
                    if key == "wave":
                        # Some records in the VSO have 'kev' which astropy
                        # doesn't recognise as a unit, so fix it.
                        waveunit = value['waveunit']
                        waveunit = 'keV' if waveunit == 'kev' else waveunit

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
        all_cols = set(data.colnames)
        first_names = [n for n in ['Start Time', 'End Time', 'Source', 'Instrument', 'Type', 'Wavelength']
                       if n in all_cols]
        extra_cols = all_cols.difference(first_names)
        all_cols = first_names + list(extra_cols)
        vsotable = data[[col for col in all_cols if data[col] is not None]]
        empty_cols = [col.info.name for col in data.itercols() if col.info.dtype.kind == 'O' and all(val is None for val in col)]
        vsotable.remove_columns(empty_cols)
        return vsotable

    def add_error(self, exception):
        self.errors.append(exception)

    def response_block_properties(self):
        """
        Returns a set of class attributes on all the response blocks.

        Returns
        -------
        s : `set`
            List of strings, containing attribute names in the response blocks.
        """
        s = {a if not a.startswith('_') else None for a in dir(self[0])}
        for resp in self[1:]:
            if len(s) == 0:
                break
            s = s.intersection({a if not a.startswith('_') else None for a in dir(resp)})

        if None in s:
            s.remove(None)
        return s
