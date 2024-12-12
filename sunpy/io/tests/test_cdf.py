
import numpy as np

import astropy.units as u

import cdflib

from sunpy.data.test import get_test_filepath
from sunpy.io._cdf import read_cdf
from sunpy.timeseries import GenericTimeSeries

filepath = get_test_filepath('solo_L2_epd-ept-north-hcad_20200713_V02.cdf')


def test_read_cdf():
    all_ts = read_cdf(filepath)
    assert isinstance(all_ts, list)
    assert len(all_ts) == 3

    ts = all_ts[0]
    print(ts.columns)
    assert isinstance(ts, GenericTimeSeries)

    col = ts.quantity('Electron_Flux_0')
    assert col.unit == u.Unit("1 / (cm2 MeV s sr)")
    # Check that fillvals are replaced by NaN
    assert np.sum(np.isnan(col)) == 189

<<<<<<< Updated upstream

def test_generic_timeseries_columns():
    # Read the CDF file using cdflib to get variable attributes
    cdf = cdflib.CDF(filepath)
    var_keys = cdf.cdf_info().rVariables + cdf.cdf_info().zVariables
    var_attrs = {key: cdf.varattsget(key) for key in var_keys}
    
    # Identify variables with 'FILLVAL'
    fillval_vars = {key: attrs['FILLVAL'] for key, attrs in var_attrs.items() if 'FILLVAL' in attrs}
    
    # Read the data using read_cdf
    all_ts = read_cdf(filepath)
    assert isinstance(all_ts, list)
    assert len(all_ts) >= 1  # Ensure at least one time series is returned
    
    # Process each GenericTimeSeries
    for ts in all_ts:
        assert isinstance(ts, GenericTimeSeries)
        
        # Map CDF variable keys to column names in GenericTimeSeries
        # Assume that column names are either the variable key or var_key_i for multi-column variables
        for col_name in ts.columns:
            # Determine the original variable key
            if '_' in col_name:
                var_key = '_'.join(col_name.split('_')[:-1])
            else:
                var_key = col_name
            
            # Skip if the variable key is not in var_attrs (可能被跳过)
            if var_key not in var_attrs:
                continue
            
            # Get the data and its attributes
            col_data = ts.quantity(col_name)
            attrs = var_attrs[var_key]
            
            # Check based on 'FILLVAL' presence
            if var_key in fillval_vars:
                # Original fill value
                fillval = fillval_vars[var_key]
                
                # Get original data from CDF
                original_data = cdf.varget(var_key)
                
                # For multi-column variables, extract the specific column
                if original_data.ndim == 2:
                    col_index = int(col_name.split('_')[-1])
                    original_col_data = original_data[:, col_index]
                else:
                    original_col_data = original_data
                
                # Positions where original data had fillval
                fill_positions = original_col_data == fillval
                
                # Ensure these positions are NaN in the column data
                assert np.all(np.isnan(col_data[fill_positions]))
                
                # Optionally, ensure no additional NaNs are introduced
                # This depends on the data; might not always be true
                # non_fill_positions = ~fill_positions
                # assert np.sum(np.isnan(col_data[non_fill_positions])) == np.sum(np.isnan(original_col_data[non_fill_positions]))
            else:
                # Ensure no NaNs are introduced unless already present
                # Check if original data had NaNs
                original_data = cdf.varget(var_key)
                
                if original_data.ndim == 2:
                    col_index = int(col_name.split('_')[-1])
                    original_col_data = original_data[:, col_index]
                else:
                    original_col_data = original_data
                
                # If original data did not have NaNs, the column data should not have NaNs
                if not np.any(np.isnan(original_col_data)):
                    assert not np.any(np.isnan(col_data))
                else:
                    # If original data had NaNs, ensure they are preserved
                    nan_positions = np.isnan(original_col_data)
                    assert np.all(np.isnan(col_data[nan_positions]))
=======
def test_check_nan_values():
    all_ts = read_cdf(filepath)
    assert isinstance(all_ts, list)

    assert len(all_ts) > 0

    for i, ts in enumerate(all_ts):
        print(f"Checking Timeseries {i + 1} for NaN values...")
        has_nan = False

        data = ts.to_dataframe()
        for column in data.columns:
            if data[column].isna().any():
                has_nan = True
                print(f"  Column '{column}' contains NaN values.")

        if not has_nan:
            print(f"  No NaN values found in Timeseries {i + 1}.")
        print()    
>>>>>>> Stashed changes
