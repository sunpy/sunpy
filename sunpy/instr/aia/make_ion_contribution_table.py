

import chianti.core as ch
import chianti.constants as c
import pandas as pd
import time


def save_contribution_csv(channel_list, temperature, density):
    """

     iterate through each column filling in the values for the contribution array with index = wavelength


    :return:
    """

    start_time = time.time()

    # define empty dataframe to fill
    data = pd.DataFrame()  # columns = ion_array, index = ion_wvl_range ,works but doesn't load != match self.wvl

    for channel in channel_list:

        # search_interval = a wavelength band around the channel center
        ion_wvl_range = []
        for i in [-2.5, 2.5]:
            ion_wvl_range.append(i + channel)
        # print(ion_wvl_range)

        # save all ions used in chiantipy
        ion_array = []
        for element in c.El:
            for level in range(1, 38, 1):
                ion_array.append(element + '_' + str(level))

        for ion in ion_array:
            try:
                # create an ion object
                ion_object = ch.ion(ion, temperature=temperature, eDensity=1.e+9)  # try without: , em=1.e+27)
                # print('ion object', ion_object)

                # find  ion intensities in the desired wavelength range
                ion_object.intensity(ion_wvl_range)

                # only do calculations for ions that have intensities in desired wavelength range
                if len(ion_object.Intensity['wvl']) != 0:
                    print(ion, len(ion_object.Intensity['wvl']))

                    # get contribution function for ion object
                    ion_object.gofnt(wvlRange=ion_wvl_range, top=1, verbose=0, plot=False)
                    contribution_function = ion_object.Gofnt['gofnt']
                    wvl = ion_object.Gofnt['wvl']

                    # round the output of contribution function to match interval wavelength range on ions
                    wavelength = round(wvl, 3)
                    # print(wavelength, ion)
                    # print(len(self.contribution_function))

                    # store ion contribution in data table based on temperature and wavelength
                    for temp, ion_contribution in zip(temperature, contribution_function):
                        # data[wavelength] = self.contribution_function
                        print(data.loc[temp, wavelength], type(data.loc[temp,wavelength]))

                        data.loc[temp, wavelength] = ion_contribution


            except (AttributeError, KeyError, UnboundLocalError, IndexError):
                # doesn't save ions if no lines found in selected interval
                #                      missing information (Elvlc file missing)
                #                      cause errors in ChiantyPy (Zion2Filename, self.Ip)
                # print('passing ', i, ) #' due to error in ChiantiPy')
                pass
                # print('ion', len(ion_contributions_array))


                # print('data:',data)

    # data in the shape G[[Temperature,wavelength of ion contribution]]
    # testing code:
    print('This code took ', time.time() - start_time, ' to run.')

    data.to_csv('all_ion_contributions.csv')



