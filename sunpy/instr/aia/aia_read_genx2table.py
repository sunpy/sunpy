from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

__author__ = "Tessa D. Wilkinson"

"""

This function reads files in a .genx data directory.
It specifically searches for instrument files and gathers information to be stored in an astropy table or dataframe.

"""

 np
import os
from scipy.io import readsav
from astropy.table import Table, Column, QTable



def aia_instr_properties_to_table(input_directory, channel_list, version, save=True):
    """

    Given  a .genx directory, this function searches for aia instrument files and saves the properties to a table.

    Parameters
    ----------
    input_directory : string,
        the path location that leads to ssw_aia_response_data/ (output of .genx.tar file)

    channel_list : list
        a list of channel wavelengths to search for

    properties: list
        a list of properties to search for

    version: int
        the version of aia instrument files to grab

    save : bool
        if True, will save a table containing all of the extracted instrument data to 'channel_properties (version #) .csv'

    Returns:
    --------

    outfile: channel_properties_[version number] .csv


    """

    check = []
    file_paths = []

    datatypes = [('wave', 'object'), ('effarea', 'object'),
                 ('geoarea', 'float32'), ('platescale', 'float32'),
                 ('wavemin', 'float32'), ('wavestep', 'float32'), ('wavenumsteps', '>i4'),
                 ('elecperdn', 'object'), ('elecperev', 'float32'),
                 ('fp_filter', 'object'), ('ent_filter', 'object'),
                 ('primary', 'object'), ('secondary', 'object'),
                 ('ccd', 'object'), ('contam', 'object')
                 ]

    properties = ['wave', 'effarea',
                  'geoarea', 'platescale',
                  'wavemin', 'wavestep', 'wavenumsteps',
                  'elecperdn', 'elecperev',
                  'fp_filter', 'ent_filter',
                  'primary', 'secondary',
                  'ccd', 'contam']

    # less ambiguous naming and manually sorted to match sorted(properties) order
    new_prop = ['quantum_efficiency_ccd', 'ccd_contamination', 'effective_area', 'electron_per_dn',
                'electron_per_ev', 'entire_filter_efficiency', 'focal_plane_filter_efficiency', 'geometric_area_ccd',
                'plate_scale', 'primary_mirror_reflectance', 'secondary_mirror_reflectance',
                'wavelength_range', 'minimum_wavelength', 'number_wavelength_intervals', 'wavelength_interval'
                ]  # efficiency or transmittance?

    # make sure data types and properties are sorted in the same order
    dt = []
    for i in sorted(properties):
        for j in datatypes:
            if j[0] == i:
                dt.append(j[1])

    table = QTable(names=new_prop, dtype=dt)

    # search through directory for instrument files for each version
    for root, dirs, files in os.walk(input_directory, topdown=False):
        for channel_file in files:
            if (os.path.join(root, channel_file)) not in check and channel_file.endswith('_fullinst') and str(
                version) in str(channel_file):
                check.append(os.path.join(root, channel_file))
                file_paths.append(os.path.join(root, channel_file))

    # search through each instrument file for channels
    for instr_file in file_paths:
        data = readsav(instr_file)['data']
        for channel_file in data.dtype.names:
            row = []

            # target number in filename that matches channels
            if channel_file.startswith('A') and channel_file.endswith('_FULL') and channel_file.find('THICK') < 0:
                start = channel_file.find('A')
                end = channel_file.find('_F')
                channel = channel_file[start + 1:end]

                if int(channel) in channel_list:
                    # match files to listed properties
                    variables = set(data[channel_file][0].dtype.fields.keys()) & set(sorted(properties))

                    # save all properties for a channel into an array
                    for n, inst_property in enumerate(sorted(variables)):
                        value = data[channel_file][0][inst_property][0]

                        row.append(value)

                    # put channel array into the table
                    table.add_row(row)

    # Make the first table column listing the channels
    C = Table.Column(data=channel_list, name='channel', dtype=int)
    table.add_column(C, index=0)
    table.add_index(['channel'])

    # could add save option!
    # if save:
    #     table.write('channel_properties_' + str(version) + '.csv')
    # return ('channel_properties_' + str(version) + '.csv')

    return table
