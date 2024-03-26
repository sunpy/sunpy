import astropy.units as u
from astropy.table import Table

from sunpy import config
from sunpy.net.base_client import BaseQueryResponse
from sunpy.time import TimeRange, parse_time
from sunpy.util.decorators import deprecated
from .table_response import iter_sort_response

__all__ = ['QueryResponse']

TIME_FORMAT = config.get("general", "time_format")

@deprecated(since="6.0", alternative="sunpy.net.vso.table_response.QueryResponse")
class QueryResponse(BaseQueryResponse):
    """
    A container for VSO Records returned from VSO Searches.
    """

    def __init__(self, lst, queryresult=None):
        super().__init__()
        self._data = lst
        self.queryresult = queryresult
        self.errors = []

        # import here to prevent circular import
        from .vso import VSOClient

        self._client = VSOClient()

    def __getitem__(self, item):
        # Always index so a list comes back
        if isinstance(item, int):
            item = slice(item, item + 1)
        return type(self)(self._data[item], queryresult=self.queryresult)

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        yield from self._data

    @property
    def blocks(self):
        return self._data

    @property
    def client(self):
        return self._client

    @client.setter
    def client(self, client):
        self._client = client

    @classmethod
    def create(cls, queryresult):
        res = list(iter_sort_response(queryresult))
        return cls(res, queryresult)

    def total_size(self):
        """ Total size of data in KB. May be less than the actual
        size because of inaccurate data providers."""
        # Warn about -1 values?
        return sum(record.size for record in self if record.size > 0)

    def time_range(self):
        """ Return total time-range all records span across. """
        return TimeRange(min(record.time.start for record in self if record.time.start is not None),
                         max(record.time.end for record in self if record.time.end is not None))

    def build_table(self):
        """
        Create a human readable table.

        Returns
        -------
        `astropy.table.QTable`
        """
        keywords = ['Start Time', 'End Time', 'Source', 'Instrument', 'Type', 'Wavelength']
        record_items = {}
        for key in keywords:
            record_items[key] = []

        def validate_time(time):
            # Handle if the time is None when coming back from VSO
            if time is None:
                return ['None']
            if record.time.start is not None:
                return [parse_time(time).strftime(TIME_FORMAT)]
            else:
                return ['N/A']

        for record in self:
            record_items['Start Time'] += validate_time(record.time.start)
            record_items['End Time'] += validate_time(record.time.end)
            record_items['Source'].append(str(record.source))
            record_items['Instrument'].append(str(record.instrument))
            if hasattr(record, 'extent') and record.extent is not None:
                record_items['Type'].append(str(record.extent.type)
                                            if record.extent.type is not None else ['N/A'])
            else:
                record_items['Type'].append('N/A')
            # If we have a start and end Wavelength, make a quantity
            if hasattr(record, 'wave') and record.wave.wavemin and record.wave.wavemax:
                unit = record.wave.waveunit
                # Convert this so astropy units parses it correctly
                if unit == "kev":
                    unit = "keV"
                record_items['Wavelength'].append(u.Quantity([float(record.wave.wavemin),
                                                              float(record.wave.wavemax)],
                                                             unit=unit))
            # If not save None
            else:
                record_items['Wavelength'].append(None)
        # If we have no wavelengths for the whole list, drop the col
        if all([a is None for a in record_items['Wavelength']]):
            record_items.pop('Wavelength')
            keywords.remove('Wavelength')
        else:
            # Make whole column a quantity
            try:
                with u.set_enabled_equivalencies(u.spectral()):
                    record_items['Wavelength'] = u.Quantity(record_items['Wavelength'])
            # If we have mixed units or some Nones just represent as strings
            except (u.UnitConversionError, TypeError):
                record_items['Wavelength'] = [str(a) for a in record_items['Wavelength']]

        return Table(record_items)[keywords]

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
