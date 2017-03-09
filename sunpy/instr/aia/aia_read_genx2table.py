"""
Reads files in a .genx data directory.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.table import QTable
import astropy.units as u

from sunpy.io.special import read_genx

__author__ = ["Tessa D. Wilkinson", "Will Barnes"]


def aia_instr_properties_to_table(channel_list, instrument_files):
    """
    Given  a .genx directory, this function searches for AIA
    instrument files and saves the properties to a table.

    Parameters
    ----------
    channel_list : array-like
        channel wavelengths to search for
    instrument_files : `list`
        AIA .genx instrument files

    Returns:
    --------
    table : `~astropy.table.QTable`
    """

    # correspondence between property names
    properties = [('wave', 'wavelength'), ('wavemin', 'minimum_wavelength'),
                  ('wavestep', 'wavelength_interval'),
                  ('wavenumsteps', 'number_wavelength_intervals'),
                  ('effarea', 'effective_area'),
                  ('geoarea', 'geometric_area_ccd'),
                  ('platescale', 'plate_scale'),
                  ('elecperdn', 'electron_per_dn'),
                  ('elecperev', 'electron_per_ev'),
                  ('fp_filter', 'focal_plane_filter_efficiency'),
                  ('ent_filter', 'entrance_filter_efficiency'),
                  ('primary', 'primary_mirror_reflectance'),
                  ('secondary', 'secondary_mirror_reflectance'),
                  ('ccd', 'quantum_efficiency_ccd'),
                  ('contam', 'ccd_contamination')]
    # corresponding units for each field
    units = [u.angstrom, u.angstrom, u.angstrom, u.dimensionless_unscaled,
             u.cm**2, u.cm**2, u.steradian/u.pixel, u.electron/u.count,
             u.electron/u.eV, u.dimensionless_unscaled,
             u.dimensionless_unscaled, u.dimensionless_unscaled,
             u.dimensionless_unscaled, u.dimensionless_unscaled,
             u.dimensionless_unscaled]
    units = {p[1]: u for p, u in zip(properties, units)}
    units['channel'] = u.angstrom
    # field name format
    field_name = 'A{0}_FULL'

    # read in values
    rows = []
    for instr_file in instrument_files:
        instrument_data = read_genx(instr_file)
        for channel in channel_list:
            if field_name.format(channel) not in instrument_data.keys():
                continue
            row = {'channel': channel}
            channel_data = instrument_data[field_name.format(channel)]
            for prop in properties:
                if prop[0].upper() not in channel_data.keys():
                    print('Cannot find {0} for channel {1} in file {2}'.format(
                        prop[0], channel, instr_file))
                    print('Setting {} to 1'.format(prop[1]))
                    row[prop[1]] = 1
                else:
                    row[prop[1]] = channel_data[prop[0].upper()]
            rows.append(row)

    # assign units
    table = QTable(rows=rows,
                   names=tuple(['channel']+[p[1] for p in properties]))
    for name in table.colnames:
        try:
            table[name].unit = units[name]
        except TypeError:
            # weird astropy table exception, units still set in all cases
            # exception seems to be only thrown when reading in the UV channels
            pass

    return table
