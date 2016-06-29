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


def aia_instr_properties_to_table(input_directory, channel_list, properties, version, save = True):
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
                    channel_information = Column(name = int(channel), data=array)

                    if int(channel) == 1600 or int(channel) == 1700:
                        pass # TODO: Fix ValueError: Inconsistent data column lengths: set([17,19])
                    #    table.add_column(channel_information)
                    else:
                        table.add_column(channel_information)

    table[properties.index('geoarea')].unit = u.cm**2

    # Last column gives indices and name of rows
    indices = Column(name='properties', data=properties)
    table.add_column(indices)
    # not sure if this works   vvv
    table.add_index(['properties'])

    assert len(table) != 0, 'Data Frame is not loading from file.'
    if save:
        table.write('channel_properties_' + str(version) + '.csv', format ='csv')

    return('channel_properties_' + str(version) + '.csv')







def aia_instr_properties_to_table_by_row(input_directory, channel_list, properties, version, save = True):
    """
    Loading in by row may be easier for obtaining data type and units.
    - WORK IN PROGRESS -
    """

    assert type(channel_list) == list
    assert type(properties) == list

    check = []
    file_paths = []

    table = Table(names = properties)

    # search through directory for instrument files for each version
    for root, dirs, files in os.walk(input_directory, topdown=False):
        for channel_file in files:
            if (os.path.join(root, channel_file)) not in check and channel_file.endswith('_fullinst') and str(version) in str(channel_file):
                # make sure not duplicated
                check.append(os.path.join(root, channel_file))
                # store in data_row
                file_paths.append(os.path.join(root, channel_file))

    assert len(file_paths) != 0, 'Did not find aia_' + str(version) + '_fullinst files in directory.'


    # search through each instrument file for channels
    for instr_file in file_paths:
        # access np.recarray from .genx file
        data = readsav(instr_file)['data']

        # pick out instrument files inside np.recarray
        for channel_file in data.dtype.names:
            row = []

            # target number in filename that matches channels
            if channel_file.startswith('A') and channel_file.endswith('_FULL') and channel_file.find('THICK') < 0:
                start = channel_file.find('A')
                end = channel_file.find('_F')
                channel = channel_file[start + 1:end]
                # get data for specific channel
                if int(channel) in channel_list:

                    # # match the properties in the rec.data_row to desired properties
                    variables = set(data[channel_file][0].dtype.fields.keys()) & set(properties)

                    # row_data
                    for inst_property in variables:

                        value = data[channel_file][0][inst_property][0]
                        if inst_property == 'contam':
                            dtype = type(value[0])
                            units = None

                            # row.append(value[0])
                            row.append(value[0])
                        else:
                            dtype = type(value)
                            if inst_property == 'geoarea':
                                units = u.cm ** 2
                            elif inst_property.startswith('wave'):
                                units = u.angstrom
                            else:
                                units = None
                            # row.append(data[value])
                            # row.append((channel, value, units,  dtype))
                            row.append(value)
                print(row)
                table.add_row(row)




    # print(data_row)
    # t = Table(rows = data_row, names = properties)
    # t = Table(rows = row)
    print(table)


    assert len(t) != 0, 'Empty Dataframe: Data Frame is not loading from file.'
    if save:
        table.write('channel_properties_' + str(version) + '.csv')
    return('channel_properties_' + str(version) + '.csv')




def aia_instr_properties_to_table_with_dict(input_directory, channel_list, properties, version, save = True):
    """
    Testing out different versions to see which is easier for obtaining data type and units.
    - WORK IN PROGRESS -
    """

    assert type(channel_list) == list
    assert type(properties) == list

    check = []
    file_paths = []

    # search through directory for instrument files for each version
    for root, dirs, files in os.walk(input_directory, topdown=False):
        for channel_file in files:
            if (os.path.join(root, channel_file)) not in check and channel_file.endswith('_fullinst') and str(version) in str(channel_file):
                # make sure not duplicated
                check.append(os.path.join(root, channel_file))
                # store in data_row
                file_paths.append(os.path.join(root, channel_file))

    assert len(file_paths) != 0, 'Did not find aia_' + str(version) + '_fullinst files in directory.'

    out_dictionary = {}

    # search through each instrument file for channels
    for instr_file in file_paths:
        # access np.recarray from .genx file
        data = readsav(instr_file)['data']

        # pick out instrument files inside np.recarray
        for channel_file in data.dtype.names:

            # target number in filename that matches channels
            if channel_file.startswith('A') and channel_file.endswith('_FULL') and channel_file.find('THICK') < 0:
                start = channel_file.find('A')
                end = channel_file.find('_F')
                channel = channel_file[start + 1:end]

                # get data for specific channel
                if int(channel) in channel_list:

                    # # match the properties in the rec.data_row to desired properties
                    variables = set(data[channel_file][0].dtype.fields.keys()) & set(properties)
                    for inst_property in variables:
                        value = data[channel_file][0][inst_property][0]
                        if inst_property == 'contam':

                            if int(channel) not in out_dictionary.keys():
                                out_dictionary[str(channel)] = [value[0]]
                            else:
                                stack = np.hstack((out_dictionary[int(channel)], [value[0]]))
                                out_dictionary[str(channel)] = stack
                            # row.append(value[0])
                        #     data_row.append((channel, inst_property, value[0], units, dtype))
                        else:

                            if int(channel) not in out_dictionary.keys():
                                out_dictionary[str(channel)] = [value]
                            else:
                                stack = np.hstack((out_dictionary[int(channel)], [value]))
                                out_dictionary[str(channel)] = stack
                        #     row.append(data[value)
                        #     data_row.append((channel, inst_property, value, units,  dtype))

    # print(row)
    # t = Table(rows = data_row, names = ('key','name', 'value', 'units', 'dtype'))
    # t = Table(rows = row)
    # print(t)
    print(out_dictionary)
    t = Table(out_dictionary, names = channel_list)


    assert len(t) != 0, 'Data Frame is not loading from file.'
    if save: # TODO: fix TypeError: unhashable type: numpy.ndarray'
        t.write('channel_properties_' + str(version) + '.csv')
    return('channel_properties_' + str(version) + '.csv')

#




def aia_instr_properties_to_dataframe(input_directory, channel_list, properties, version, save = True):
    """
    This definition reads the instrument file aia_V6_all_fullinst, which was obtained from ssw_aia_response_data inside
    ssw_aia_response_genx.tar.gz (an output file saved from SolarSoft Ware (SSW)). It will extract the instrument data and save
    it a dataframe file for easier access.
    Returns
    -------
    :returns outfile to dataframe
        Each column lists the channel information and each row is a different peroperty from aia_inst_genx.
    :returns 'channel_properties.csv', dataframe with wavelength centers for each channel is a column and listed in the column are the properties from the is returned where it has keys of the properties of the channel
    Notes:
    np.recarray store information with shape(1,0) and are quite nested.
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
                # store in data_row
                file_paths.append(os.path.join(root, name))

    assert len(file_paths) != 0, 'Did not find aia_' + str(version) + '_fullinst files in directory.'

    # store in dataframe
    df = pd.DataFrame()

    for instr_file in file_paths:
        # access np.recarray from .genx file
        ssw_array = readsav(instr_file)
        data = ssw_array['data']

        # pick out instrument files inside np.recarray
        for name in data.dtype.names:

            # Read in data from Channels from  A##_FULL  which match format in UV and EUV files
            if name.startswith('A') and name.endswith('_FULL') and name.find('THICK') < 0:

                # target number in filename that matches wavelength_center
                start = name.find('A')
                end = name.find('_F')
                wavelength = name[start + 1:end]

                # UV files have 'thick' files that need to be filtered out for no duplicates
                if int(wavelength) in channel_list:

                    # match the properties in the rec.array to desired properties
                    variables = set(data[name][0].dtype.fields.keys()) & set(properties)

                    # then iterate through the properties to fill a dataframe with rec.array information
                    for property in variables:
                        try:
                            df.loc[property, wavelength] = data[name][0][property][0]
                        # contam and crossarea are sequences -- # TODO: work on this try statement
                        except ValueError:
                            df.loc[property, wavelength] = data[name][0][property]
                else:
                    pass

    assert len(df) != 0, 'Data Frame is not loading from file.'

    if save:
        #save to .csv outfile
        df.to_csv('channel_properties_' + str(version) + '.csv')

    # WANT:
    # return ('channel_properties_' + str(version) + '.csv')
    return df

# I can use this dataframe and load in to table with chosen dtypes and units. not generic  though
