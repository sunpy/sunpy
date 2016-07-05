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
from astropy.table import Table, Column
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

    outfile: channel_properties[version number] .csv

    Notes;

    astropy Table was a recommneded data structre to house this content.
    I am not as familiar with it as compared to dataframes.
    It stores items by column, which is good for us to store things by channel.
    I am not as familiar as to how to apply the dtype, as it doens't like the rec.array dtype format.
    I also am unsure of how to load just one index from this table...

    """

    assert type(channel_list) == list
    assert type(properties) == list

    check = []
    file_paths = []

    # search through directory for instrument files for each version
    for root, dirs, files in os.walk(input_directory, topdown=False):
        for name in files:
            if (os.path.join(root, name)) not in check and name.endswith('_fullinst') and str(version) in str(name):
                # make sure not duplicated
                check.append(os.path.join(root, name))
                # store in array
                file_paths.append(os.path.join(root, name))

    assert len(file_paths) != 0, 'Did not find aia_' + str(version) + '_fullinst files in directory.'

    # # store in table
    table = Table()

    # search through each instrument file for channels
    for instr_file in file_paths:
        # access np.recarray from .genx file
        data = readsav(instr_file)['data']

        # pick out instrument files inside np.recarray
        for name in data.dtype.names:
            array = []

            # target number in filename that matches channels
            if name.startswith('A') and name.endswith('_FULL') and name.find('THICK') < 0:
                start = name.find('A')
                end = name.find('_F')
                channel = name[start + 1:end]
                # get data for specific channel
                if int(channel) in channel_list:
                    # # match the properties in the rec.array to desired properties
                    variables = set(data[name][0].dtype.fields.keys()) & set(properties)

                    # iterate through the properties and append to array
                    for inst_property in variables:
                        # NOTE: ValueError: setting an array element with a sequence solved by making each element in array a string but I loose type information...
                        if inst_property == 'contam':
                            array.append(str(data[name][0][inst_property][0][0]))
                        else:
                            array.append(str(data[name][0][inst_property][0]))

                    # dtype = []
                    # for i in array:
                    #     print(i, type(i))
                    #     dtype.append(type(i))
                    # create column of property information per channel
                    # TODO: implement dtype, unit
                    channel_information = Column(name=int(channel), data=array)

                    if int(channel) == 1600 or int(channel) == 1700:
                        pass  # TODO: Fix ValueError: Inconsistent data column lengths: set([17,19])
                    # table.add_column(channel_information)
                    else:
                        table.add_column(channel_information)

    # Last column gives indices and name of rows
    indices = Column(name='properties', data=properties)
    table.add_column(indices)
    # not sure if this works   vvv
    table.add_index('properties')
    # print(table)
    assert len(table) != 0, 'Data Frame is not loading from file.'
    if save:
        table.write('channel_properties_' + str(version) + '.csv', format='csv')

    return ('channel_properties_' + str(version) + '.csv')
