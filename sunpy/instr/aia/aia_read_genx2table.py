from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

__author__ = "Tessa D. Wilkinson"

"""

This function reads files in a .genx data directory.
It specifically searches for instrument files and gathers information to be stored in an astropy table or dataframe.

"""

import numpy as np
import os
from scipy.io import readsav
from astropy.table import Table, Column, QTable
import astropy.units as u
import pandas as pd


def aia_instr_properties_to_table(input_directory, channel_list, properties, version, save=True):
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

    Notes;


    """

    assert type(channel_list) == list
    assert type(properties) == list

    check = []
    file_paths = []

    # table = Table(names = properties)
    # table = Table(names=properties, dtype = ('int32', 'float32', 'int16', 'float32', '>f4', 'int16', '>f4', 'float32', '>f4', '>f4', 'float32', 'int16', '>f4', 'int16', '>f4' ,'int16' ,'float32', 'float32', '>f4', '>f4'))

    # used this one!
    #datatypes = [('wavenumsteps','>i4'), ('geoarea','float32'),('usecontam', 'int16'),('elecperdn',  '>f4'), ('contam','>f4'), ('numfilters','int16'), ('primary', '>f4'), ('wavemin', 'float32'), ('wave','>f4'), ('ent_filter', '>f4'), ('wavestep', 'float32'), ('useerror','int16'), ('effarea', '>f4'),('wavelog', 'int16'),('ccd', '>f4') ,('contamthick', 'int16') ,('platescale', 'float32'), ('elecperev','float32'),('fp_filter', '>f4')]  #, ('secondary','>f4'), ('elecperphot', 'float32')
    datatypes = [('wavenumsteps','>i4'), ('geoarea','float32'),('usecontam', 'int16'),('elecperdn',  'object'), ('contam','object'), ('numfilters','int16'), ('primary', 'object'), ('wavemin', 'float32'), ('wave','object'), ('ent_filter', 'object'), ('wavestep', 'float32'), ('useerror','int16'), ('effarea', 'object'),('wavelog', 'int16'),('ccd', 'object') ,('contamthick', 'int16') ,('platescale', 'float32'), ('elecperev','float32'),('fp_filter', 'object'), ('secondary','object')]#, ('elecperphot', 'float32')]

    # datatypes = ('int32', 'float32', 'int16', 'float32', '>f4', 'int16', '>f4', 'float32', '>f4', '>f4', 'float32', 'int16', '>f4', 'int16', '>f4' ,'int16' ,'float32', 'float32',  '>f4', '>f4')

    # print(len(properties))
    # print(len(datatypes))

    dt = []
    for i in sorted(properties):
        for j in datatypes:
            if j[0] == i:
                # print(i,j)
                dt.append(j[1])

    table = QTable(names = sorted(properties), dtype = dt)


    # search through directory for instrument files for each version
    for root, dirs, files in os.walk(input_directory, topdown=False):
        for channel_file in files:
            if (os.path.join(root, channel_file)) not in check and channel_file.endswith('_fullinst') and str(
                version) in str(channel_file):
                check.append(os.path.join(root, channel_file))
                file_paths.append(os.path.join(root, channel_file))

    assert len(file_paths) != 0, 'Did not find aia_' + str(version) + '_fullinst files in directory.'


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
                        dtype = value.dtype
                        if str(dtype) == '>f4':
                            row.append(value)
                            # print(n, inst_property, dtype)
                        else:
                            row.append(value)
                            # print(n, inst_property, dtype)

                    # put channel array into the table
                    table.add_row(row)


    # Make the first table column listing the channels
    C = Table.Column(data=channel_list, name='channel', dtype=int)
    table.add_column(C, index = 0)
    table.add_index(['channel'])

    # print('OUT- TABLE------------------', table)

    assert len(table) != 0, 'Empty Dataframe: Data Frame is not loading from file.'
    # if save:
    #     table.write('channel_properties_' + str(version) + '.csv')
    # return ('channel_properties_' + str(version) + '.csv')

    #^^^^^^^^^Returned TypeError: unhashable type: 'numpy.ndarray'
    # so maybe it is better to read in each time

    return table
