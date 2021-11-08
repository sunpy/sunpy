import cdflib
import pandas as pd
from cdflib.epochs import CDFepoch

import astropy.units as u

from sunpy import log
from sunpy.timeseries import GenericTimeSeries

__all__ = ['read_cdf']


def read_cdf(fname):
    """
    Read a CDF file that follows the ISTP/IACG guidelines.

    Parameters
    ----------
    fname : path-like
        Location of single CDF file to read.

    Returns
    -------
    list[GenericTimeSeries]
        A list of time series objects, one for each unique time index within
        the CDF file.

    References
    ----------
    Space Physics Guidelines for CDF https://spdf.gsfc.nasa.gov/sp_use_of_cdf.html
    """
    cdf = cdflib.CDF(str(fname))

    # Extract the time varying variables
    cdf_info = cdf.cdf_info()
    meta = cdf.globalattsget()
    all_var_keys = cdf_info['rVariables'] + cdf_info['zVariables']
    var_attrs = {key: cdf.varattsget(key) for key in all_var_keys}
    # Get keys that depend on time
    var_keys = [var for var in var_attrs if 'DEPEND_0' in var_attrs[var]]

    # Get unique time index keys
    time_index_keys = sorted(set([var_attrs[var]['DEPEND_0'] for var in var_keys]))

    all_ts = []
    # For each time index, construct a GenericTimeSeries
    for index_key in time_index_keys:
        index = cdf.varget(index_key)
        # TODO: use to_astropy_time() instead here when we drop pandas in timeseries
        index = CDFepoch.to_datetime(index)
        df = pd.DataFrame(index=pd.DatetimeIndex(name=index_key, data=index))
        units = {}

        for var_key in sorted(var_keys):
            attrs = var_attrs[var_key]
            if attrs['DEPEND_0'] != index_key:
                continue

            # Get data
            if cdf.varinq(var_key)['Last_Rec'] == -1:
                log.debug(f'Skipping {var_key} in {fname} as it has zero elements')
                continue

            data = cdf.varget(var_key)
            # Get units
            unit = attrs['UNITS']
            if unit in ['None', '', 'unitless']:
                unit = u.dimensionless_unscaled
            else:
                unit = u.Unit(unit)

            if data.ndim == 2:
                # Multiple columns, give each column a unique label
                for i, col in enumerate(data.T):
                    df[var_key + f'_{i}'] = col
                    units[var_key + f'_{i}'] = unit
            else:
                # Single column
                df[var_key] = data
                units[var_key] = unit

        all_ts.append(GenericTimeSeries(data=df, units=units, meta=meta))

    return all_ts
