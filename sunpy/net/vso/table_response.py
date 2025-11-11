"""
Classes and helper functions for VSO responses.
"""
from collections import defaultdict
from collections.abc import Mapping

import numpy as np
from zeep.helpers import serialize_object

import astropy.units as u
from astropy.table import TableAttribute
from astropy.time import Time

from sunpy.net.base_client import QueryResponseTable
from sunpy.time import parse_time

__all__ = ['VSOQueryResponseTable']


class VSOQueryResponseTable(QueryResponseTable):
    hide_keys = ['fileid', 'fileurl', 'Info Required']
    errors = TableAttribute(default=[])
    size_column = 'Size'

    @classmethod
    def from_zeep_response(cls, response, *, client):
        """
        Construct a table response from the zeep response.
        """
        data = []
        for record in response:
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
                        # doesn't recognize as a unit, so fix it.
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
                        if key == "time" and subvalue is not None:
                            key_template = f"{subkey.capitalize()} {key.capitalize()}"
                        else:
                            key_template = f"{key.capitalize()} {subkey.capitalize()}"
                        row[key_template] = subvalue
            data.append(row)
        data = cls(data, client=client)
        # Parse times
        for col in data.colnames:
            if col.endswith('Time'):
                try:
                    # Try to use a vectorised call to parse_time
                    data[col] = parse_time(data[col])
                except Exception:
                    # If that fails, parse dates one by one. This is needed if
                    # VSO returns a variety of different date format strings
                    times = []
                    mask = []
                    for i, t in enumerate(data[col]):
                        if t is not None:
                            times.append(parse_time(t))
                        else:
                            # Create a dummy time and mask it later
                            times.append(Time(val=0, format='mjd'))
                            mask.append(i)
                    data[col] = Time(times)
                    data[col][mask] = np.ma.masked

                data[col].format = 'iso'

        # Reorder the columns to put the most useful ones first.
        return data._reorder_columns(['Start Time', 'End Time', 'Source',
                                      'Instrument', 'Type', 'Wavelength'],
                                     remove_empty=True)

    def add_error(self, exception):
        self.errors.append(exception)
