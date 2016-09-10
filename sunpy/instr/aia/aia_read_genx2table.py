"""
Reads files in a .genx data directory.
"""
from __future__ import (absolute_import, division, print_function,unicode_literals)
import os

import numpy as np
from scipy.io import readsav
from astropy.table import QTable
import astropy.units as u

__author__ = "Tessa D. Wilkinson"


def aia_instr_properties_to_table(input_directory, channel_list, version, save=True):
    """

    Given  a .genx directory, this function searches for aia instrument files and saves the properties to a table.

    Parameters
    ----------
    input_directory : string,
        the path location that leads to ssw_aia_response_data/ (output of .genx.tar file)
    channel_list : list
        a list of channel wavelengths to search for
    version: int
        the version of aia instrument files to grab

    Returns:
    --------
    outfile: channel_properties_[version number] .csv
    """

    #correspondence between property names
    properties = [('wave','wavelength'),('wavemin','minimum_wavelength'),('wavestep','wavelength_interval'),
                  ('wavenumsteps','number_wavelength_intervals'),('effarea','effective_area'),
                  ('geoarea','geometric_area_ccd'),('platescale','plate_scale'),('elecperdn','electron_per_dn'),
                  ('elecperev','electron_per_ev'),('fp_filter','focal_plane_filter_efficiency'),
                  ('ent_filter','entire_filter_efficiency'),('primary','primary_mirror_reflectance'),
                  ('secondary','secondary_mirror_reflectance'),('ccd','quantum_efficiency_ccd'),
                  ('contam','ccd_contamination')]
    #corresponding units for each field
    unitless = u.m/u.m
    units = [u.angstrom,u.angstrom,u.angstrom,unitless,u.cm**2,u.cm**2,1/u.cm,u.electron/u.count,u.electron/u.eV,
             unitless,unitless,unitless,unitless,unitless,unitless]
    units = {p[1] : u for p,u in zip(properties,units)}
    units['channel'] = u.angstrom
    #set instrument files; these are the names used by SSW, should be static
    instrument_files = [os.path.join(input_directory,'aia_V{0}_all_fullinst'.format(version)),
                        os.path.join(input_directory,'aia_V{0}_fuv_fullinst'.format(version))]
    field_name = 'A{0}_FULL'

    #read in values
    rows = []
    for instr_file in instrument_files:
        instrument_data = readsav(instr_file)['data']
        for channel in channel_list:
            if field_name.format(channel).upper() not in instrument_data.dtype.names:
                continue
            row = {'channel':channel}
            channel_data = instrument_data[field_name.format(channel)][0]
            for prop in properties:
                if prop[0].upper() not in channel_data.dtype.names:
                    print('Cannot find {0} for channel {1} in file {2}'.format(prop[0],channel,instr_file))
                    print('Setting {} to 1'.format(prop[1]))
                    row[prop[1]] = 1
                else:
                    row[prop[1]] = channel_data[prop[0]][0]
            rows.append(row)

    #assign units
    table = QTable(rows=rows,names=tuple(['channel']+[p[1] for p in properties]))
    for name in table.colnames:
        try:
            table[name].unit = units[name]
        except TypeError:
            #weird astropy table exception, units still set in all cases
            #exception seems to be only thrown when reading in the UV channels
            pass

    return table
