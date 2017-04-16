import ChiantiPy.core as ch
import ChiantiPy.tools.constants as c
import pandas as pd
import time
import numpy as np
import math


def magnitude(x):
    return int(math.log10(x))


def save_contribution_csv(channel_list, wavelength_response, temperature, density):
    """

     iterate through each column filling in the values for the contribution array with index = wavelength


    :return: a csv file
        The csv file contains the contribution function for each ion calculated in ChiantiPy. The columns are
        the wavelengths and each row is the temperature that the transition occurs.
    """

    start_time = time.time()

    # define empty dataframe to fill
    data = pd.DataFrame()  # columns = ion_array, index = ion_wvl_range ,works but doesn't load != match self.wvl

    # save all ions names used in chiantipy as strings in an array
    ion_array = []
    for element in c.El:
        for level in range(1, 38, 1):
            ion_array.append(element + '_' + str(level))

    # define wavelength range for each channel based on wavelength response values
    wavelength_range_dict = {}
    non_zero_array = []
    for wavelength_center, wr_array in wavelength_response.items():

        # return array of indices for each channel wavelength response array that are nonzero
        max = (wr_array['response'].value).max()
        if max > 0.1:
            index_array = np.where(wr_array['response'].value > 0.05)
        else:
            index_array = np.where(wr_array['response'].value > 0.015)  # smaller range than 0.01

        # define start and end indices
        index_start = index_array[0][0]
        index_end = index_array[0][-1]

        # find wavelength values for those indices
        wavelength_start = wr_array['wavelength'][index_start]
        wavelength_end = wr_array['wavelength'][index_end]
        # print(wavelength_center, wavelength_end- wavelength_start)

        # save wavelength ranges for each channel
        wavelength_range_dict[wavelength_center] = [wavelength_start, wavelength_end]

    # define density (next step to make this density = [pressure = 10 ** 15 #* (u.Kelvin / u.cm ** 3) /temperature])
    density = 1.e+9

    # iterate through each ion for each small wavelength range
    # could be more efficient with a way in chiantipy that allows for looking up transitions of ions per wavelength without making an ion object
    for ion in ion_array:
        for wavelength_center, wavelength_range in wavelength_range_dict.items():

            # define indices to probe
            ion_wvl_range = wavelength_range

            try:

                # define ion object
                ion_object = ch.ion(ion, temperature=temperature, eDensity=density)

                # find  ion intensities in the desired wavelength range
                ion_object.intensity(ion_wvl_range)

                # limit to transitions that happen in the wavelength range
                if len(ion_object.Intensity['wvl']) == 0:
                    continue

                # define max intensity for reference
                max_intensity = ion_object.Intensity['intensity'].max()
                mag_max = magnitude(max_intensity)
                # print(mag_max)

                # limit to ion transitions above threshold: 10**-28
                if mag_max <= -28:
                    continue

                # calculate contribution function
                ion_object.gofnt(wvlRange=ion_wvl_range, top=1, verbose=0, plot=False)
                contribution_function = ion_object.Gofnt['gofnt']
                # print(contribution_function)
                wvl = ion_object.Gofnt['wvl']

                # round the output wavelengths of contribution function
                wavelength = round(wvl, 3)
                # print(wavelength, ion)

                # store ion contribution in data table based on temperature and wavelength
                for temp, ion_contribution in zip(temperature, contribution_function):
                    # print(temp, ion_contribution)
                    data.loc[temp, wavelength] = ion_contribution

            except (AttributeError, KeyError, UnboundLocalError, IndexError):
                # doesn't save ions if no lines found in selected interval - ChiantiPy auto-prints message
                #                      missing information (Elvlc file missing)- ChiantiPy auto-prints message
                #                      cause errors in ChiantyPy (Zion2Filename, self.Ip)
                continue

    # testing code:
    print('This code took ', time.time() - start_time, ' seconds to run.')

    data.to_csv('test.csv')
