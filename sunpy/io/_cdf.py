"""
This module provides a CDF file reader for internal use.
"""
import inspect

import cdflib
import numpy as np
import pandas as pd
from cdflib.epochs import CDFepoch
from packaging.version import Version

import astropy.units as u

from sunpy import log
from sunpy.timeseries import GenericTimeSeries
from sunpy.util.exceptions import warn_user

__all__ = ['read_cdf']


def read_cdf(fname, **kwargs):
    """
    Read a CDF file that follows the ISTP/IACG guidelines.

    Parameters
    ----------
    fname : path-like
        Location of single CDF file to read.
    **kwargs : dict
        Additional keyword arguments are handed to ``cdflib.CDF`` reader.

    Returns
    -------
    list[GenericTimeSeries]
        A list of time series objects, one for each unique time index within
        the CDF file.

    References
    ----------
    Space Physics Guidelines for CDF https://spdf.gsfc.nasa.gov/sp_use_of_cdf.html
    """
    # Limit to kwargs that exist in cdflib.CDF
    filtered_kwargs = {k: v for k, v in kwargs.items() if k in [param.name for param in inspect.signature(cdflib.CDF).parameters.values()]}
    cdf = cdflib.CDF(str(fname), **filtered_kwargs)
    # Extract the time varying variables
    cdf_info = cdf.cdf_info()
    meta = cdf.globalattsget()
    if hasattr(cdflib, "__version__") and Version(cdflib.__version__) >= Version("1.0.0"):
        all_var_keys = cdf_info.rVariables + cdf_info.zVariables
    else:
        all_var_keys = cdf_info['rVariables'] + cdf_info['zVariables']
    var_attrs = {key: cdf.varattsget(key) for key in all_var_keys}
    # Get keys that depend on time
    var_keys = [var for var in var_attrs if 'DEPEND_0' in var_attrs[var] and var_attrs[var]['DEPEND_0'] is not None]

    # Get unique time index keys
    time_index_keys = sorted(set([var_attrs[var]['DEPEND_0'] for var in var_keys]))

    all_ts = []
    # For each time index, construct a GenericTimeSeries
    for index_key in time_index_keys:
        try:
            index = cdf.varget(index_key)
        except ValueError:
            # Empty index for cdflib >= 0.3.20
            continue
        # TODO: use to_astropy_time() instead here when we drop pandas in timeseries
        index = CDFepoch.to_datetime(index)
        df_dict = {}
        units = {}

        for var_key in sorted(var_keys):
            attrs = var_attrs[var_key]
            # If this variable doesn't depend on this index, continue
            if attrs['DEPEND_0'] != index_key:
                continue

            # Get data
            if hasattr(cdflib, "__version__") and Version(cdflib.__version__) >= Version("1.0.0"):
                var_last_rec = cdf.varinq(var_key).Last_Rec
            else:
                var_last_rec = cdf.varinq(var_key)['Last_Rec']
            if var_last_rec == -1:
                log.debug(f'Skipping {var_key} in {fname} as it has zero elements')
                continue

            data = cdf.varget(var_key)

            # Set fillval values to NaN
            # It would be nice to properley mask these values to work with
            # non-floating point (ie. int) dtypes, but this is not possible with pandas
            if np.issubdtype(data.dtype, np.floating):
                if 'FILLVAL' in attrs:
                    data[data == attrs['FILLVAL']] = np.nan


            # Get units
            if 'UNITS' in attrs:
                unit_str = attrs['UNITS']
                try:
                    unit = u.Unit(unit_str)
                except ValueError:
                    if unit_str in _known_units:
                        unit = _known_units[unit_str]
                    else:
                        warn_user(f'astropy did not recognize units of "{unit_str}". '
                                  'Assigning dimensionless units. '
                                  'If you think this unit should not be dimensionless, '
                                  'please raise an issue at https://github.com/sunpy/sunpy/issues')
                        unit = u.dimensionless_unscaled
            else:
                warn_user(f'No units provided for variable "{var_key}". '
                          'Assigning dimensionless units.')
                unit = u.dimensionless_unscaled

            if data.ndim > 2:
                # Skip data with dimensions >= 3 and give user warning
                warn_user(f'The variable "{var_key}" has been skipped because it has more than 2 dimensions, which is unsupported.')
            elif data.ndim == 2:
                # Multiple columns, give each column a unique label
                for i, col in enumerate(data.T):
                    df_dict[var_key + f'_{i}'] = col
                    units[var_key + f'_{i}'] = unit
            else:
                # Single column
                df_dict[var_key] = data
                units[var_key] = unit
        df = pd.DataFrame(df_dict, index=pd.DatetimeIndex(name=index_key, data=index))
        all_ts.append(GenericTimeSeries(data=df, units=units, meta=meta))

    if not len(all_ts):
        log.debug(f'No data found in file {fname}')
    return all_ts


# Unfortunately (unlike e.g. FITS), there is no standard for the strings that
# CDF files use to represent units. To allow for this we maintain a dictionary
# mapping unit strings to their astropy unit equivalents.
#
# Please only add new entries if
#   1. A user identifies which specific mission/data source they are needed for
#   2. The mapping from the string to unit is un-ambiguous. If we get this
#      wrong then users will silently have the wrong units in their data!
_known_units = {'ratio': u.dimensionless_unscaled,
                'NOTEXIST': u.dimensionless_unscaled,
                'Unitless': u.dimensionless_unscaled,
                'unitless': u.dimensionless_unscaled,
                'Quality_Flag': u.dimensionless_unscaled,
                'None': u.dimensionless_unscaled,
                'none': u.dimensionless_unscaled,
                ' none': u.dimensionless_unscaled,
                'counts': u.dimensionless_unscaled,
                'cnts': u.dimensionless_unscaled,

                'microW m^-2': u.mW * u.m**-2,

                'years': u.yr,
                'days': u.d,

                '#/cc': u.cm**-3,
                '#/cm^3': u.cm**-3,
                'cm^{-3}': u.cm**-3,
                'particles cm^-3': u.cm**-3,
                'n/cc (from moments)': u.cm**-3,
                'n/cc (from fits)': u.cm**-3,
                'Per cc': u.cm**-3,
                '#/cm3': u.cm**-3,
                'n/cc': u.cm**-3,

                'km/sec': u.km / u.s,
                'km/sec (from fits)': u.km / u.s,
                'km/sec (from moments)': u.km / u.s,
                'Km/s': u.km / u.s,

                'Volts': u.V,

                'earth radii': u.earthRad,
                'Re': u.earthRad,
                'Earth Radii': u.earthRad,
                'Re (1min)': u.earthRad,
                'Re (1hr)': u.earthRad,

                'Degrees': u.deg,
                'degrees': u.deg,
                'Deg': u.deg,
                'deg (from fits)': u.deg,
                'deg (from moments)': u.deg,
                'deg (>200)': u.deg,

                'Deg K': u.K,
                'deg_K': u.K,
                '#/{cc*(cm/s)^3}': (u.cm**3 * (u.cm / u.s)**3)**-1,
                'sec': u.s,
                'Samples/s': 1 / u.s,

                'seconds': u.s,
                'nT GSE': u.nT,
                'nT GSM': u.nT,
                'nT DSL': u.nT,
                'nT SSL': u.nT,
                'nT (1min)': u.nT,
                'nT (3sec)': u.nT,
                'nT (1hr)': u.nT,
                'nT (>200)': u.nT,

                'msec': u.ms,
                'milliseconds': u.ms,

                '#/cm2-ster-eV-sec': 1 / (u.cm**2 * u.sr * u.eV * u.s),
                '#/(cm^2*s*sr*MeV/nuc)': 1 / (u.cm**2 * u.s * u.sr * u.MeV),
                '#/(cm^2*s*sr*Mev/nuc)': 1 / (u.cm**2 * u.s * u.sr * u.MeV),
                '#/(cm^2*s*sr*Mev/nucleon)': 1 / (u.cm**2 * u.s * u.sr * u.MeV),
                '#/(cm2-steradian-second-MeV/nucleon) ': 1 / (u.cm**2 * u.s * u.sr * u.MeV),
                '1/(cm2 Sr sec MeV/nucleon)': 1 / (u.cm**2 * u.sr * u.s * u.MeV),
                '1/(cm**2-s-sr-MeV)': 1 / (u.cm**2 * u.s * u.sr * u.MeV),
                '1/(cm**2-s-sr-MeV/nuc.)': 1 / (u.cm**2 * u.s * u.sr * u.MeV),
                '(cm^2 s sr MeV/n)^-1': 1 / (u.cm**2 * u.s * u.sr * u.MeV),
                'cm!E-2!Nsr!E-1!Nsec!E-1!N(MeV/nuc)!E-1!N': 1 / (u.cm**2 * u.s * u.sr * u.MeV),
                'cm!E-2!Nsr!E-1!Nsec!E-1!NMeV!E-1!N': 1 / (u.cm**2 * u.s * u.sr * u.MeV),
                '1/(cm^2 sec ster MeV)': 1 / (u.cm**2 * u.s * u.sr * u.MeV),
                '(cm^2 s sr MeV)^-1': 1 / (u.cm**2 * u.s * u.sr * u.MeV),
                'cnts/sec/sr/cm^2/MeV': 1 / (u.cm**2 * u.s * u.sr * u.MeV),

                'particles / (s cm^2 sr MeV)': 1 / (u.cm**2 * u.s * u.sr * u.MeV),
                'particles / (s cm^2 sr MeV/n)': 1 / (u.cm**2 * u.s * u.sr * u.MeV),
                'particles/(s cm2 sr MeV/n)': 1 / (u.cm**2 * u.s * u.sr * u.MeV),

                '1/(cm**2-s-sr)': 1 / (u.cm**2 * u.s * u.sr),
                '1/(SQcm-ster-s)': 1 / (u.cm**2 * u.s * u.sr),
                '1/(SQcm-ster-s)..': 1 / (u.cm**2 * u.s * u.sr),

                'photons cm^-2 s^-1': 1 / (u.cm**2 * u.s),

                'Counts/256sec': 1 / (256 * u.s),
                'Counts/hour': 1 / u.hr,
                'counts/min': 1 / u.min,
                'counts / s': 1/u.s,
                'counts/s': 1/u.s,
                'cnts/sec': 1/u.s,
                'counts s!E-1!N': 1/u.s,
                }
