from __future__ import absolute_import

import tempfile

import pytest
import numpy as np

import sunpy.map
import sunpy.data.test as test


# Define the original and prepped images first so they're available to all functions

@pytest.fixture
def original():
    return sunpy.map.Map(test.get_test_filepath("aia_171_level1.fits"))

@pytest.fixture
def prep_map(original):
    return aiaprep(original)


def test_aiaprep(original, prep_map):
    # Test that header info for the map has been correctly updated
    # Check all of these for Map attributes and .meta values?
    # Check array shape
    assert prep_map.data.shape == original.data.shape
    # Check crpix values
    assert prep_map.meta['crpix1'] == prep_map.data.shape[1]/2.0 + 0.5
    assert prep_map.meta['crpix2'] == prep_map.data.shape[0]/2.0 + 0.5
    # Check cdelt values
    assert prep_map.meta['cdelt1']/0.6 == int(prep_map.meta['cdelt1']/0.6)
    assert prep_map.meta['cdelt2']/0.6 == int(prep_map.meta['cdelt2']/0.6)
    # Check rotation value, I am assuming that the inaccuracy in
    # the CROTA -> PCi_j matrix is causing the inaccuracy here
    np.testing.assert_allclose(prep_map.rotation_matrix, np.identity(2), rtol=1e-5, atol=1e-8)
    # Check level number
    assert prep_map.meta['lvl_num'] == 1.5


def test_filesave(prep_map):
    # Test that adjusted header values are still correct after saving the map
    # and reloading it.
    afilename = tempfile.NamedTemporaryFile(suffix='fits').name
    prep_map.save(afilename, clobber=True)
    load_map = sunpy.map.Map(afilename)
    # Check crpix values
    assert load_map.meta['crpix1'] == prep_map.data.shape[1]/2.0 + 0.5
    assert load_map.meta['crpix2'] == prep_map.data.shape[0]/2.0 + 0.5
    # Check cdelt values
    assert load_map.meta['cdelt1']/0.6 == int(load_map.meta['cdelt1']/0.6)
    assert load_map.meta['cdelt2']/0.6 == int(load_map.meta['cdelt2']/0.6)
    # Check rotation value
    np.testing.assert_allclose(prep_map.rotation_matrix, np.identity(2), rtol=1e-5, atol=1e-8)
    # Check level number
    assert load_map.meta['lvl_num'] == 1.5




# def test_aia_read_genx2table():
#     # Check that keywords are lists
#     assert type(channel_list) == list
#     #assert type(properties) == list   # no longer keyword - but should it be?
#     assert len(properties) == len(datatypes) == len(new_prop) # - needed?
#     #
#     # should match length
#     # print(len(datatypes))
#     # print(len(properties))
#     # print(len(new_prop))
#
#     # check that it finds files to look in the genx directory
#     assert len(file_paths) != 0, 'Did not find aia_' + str(version) + '_fullinst files in directory.'
#
#     # Check that the table is filled
#     assert len(table) != 0, 'Empty Table: Data is not loading from file.'




# def test_response():
    # Test for response function

    # Check initial properties are loaded in as arrays
    # assert type(self.channel_list) == array

    # test get wavelength range
    # if no channel specified, wavelength range should match var['wavelength_range']
    # assert here..?


    # test effective area

    # # This test confirms that I need to index the values differently than just the index of the channel
    # for value in var['fp_filter']:
    #     if 0.348 < values < 0.349: # trying to replicate table data
    #         index = var['fp_filter'].index(value) # 684   != 94
    #         print(index, value)
    #         print(self.get_wavelength_range(channel)[index]) # 93.9 !!

    # assert statement, these two should be the same
    # print(len(var['fp_filter']))
    # print(var['wavenumsteps'])

    # print what's being called, perform test if not close to what's expected
    # print(var['primary'][index], var['secondary'][index], var['fp_filter'][index], var['ent_filter'][index])
    # print(reflectance, transmission_efficiency)


    # test get_wavelength_response
    # use genx files
        # print(gain)# 18.3 for all channels from .genx file
    # after else
        # print('elecperphot: ', ev_per_wavelength)
        # ^^^ assert for wavelength 94, expect about 131.89 nm
        # these two should NOT be the same,
        # print('wave', var['wave'][index] , var['wave'][channel])

        # print('e', electron_per_ev, calc_electron_per_ev)            # these are approximately the same

        # print('epl', electron_per_wavelength)
        # print(ccdgain)

        # print('channel:', channel )
        # print('index:', index )
        # print('value:', value )
        # print('platescale: ', var['platescale'])
        # print('gain: ', [0.664, 1.785, 3.365, 1.217, 1.14, 0.041, 0.027, 0.0024, 0.0046][
        #     self.channel_list.index(channel)], calculated_gain * electron_per_dn)  # ~ 17 or 18.3 dn/phot


        # test gain:
        # # move to aia_test.pu as a reference
        # if use_table2_gain:
        #     # Uses Gain from table Table 12 to calculate instrument response
        #     # gain values from Boerner paper
        #     gain_table2 = {94: 2.128, 131: 1.523, 171: 1.168, 193: 1.024, 211: 0.946, 304: 0.658, 335: 0.596,
        #                    1600: 0.125, 1700: 0.118}
        #     self.system_gain = gain_table2[channel]


        # test calculate wavelength response:
        # if use_table2_gain:
        #     self.calculate_effective_area(channel)
        #     self.calculate_system_gain(channel, use_table2_gain=True)
        #     self.wavelength_response = self.effective_area * self.system_gain * dn_per_photon

        # test calculate wavelength response
        # elif use_response_table2:
        #     # returns the instrument response values listed in Table 2
        #     i = self.channel_list.index(channel)
        #     response = [0.664, 1.785, 3.365, 1.217, 1.14, 0.041, 0.027, 0.0024, 0.0046][i]
        #     self.wavelength_response = response

        # temperature response:

        # print('matrix:' , len(matrix)) # 61
        # print(discrete_wavelengths) # returns array   (1,3) for channel 94, only three ions probed - check

        # check that these are the same length -- good! 8751 == 8751
        # print( len(self.wavelength_response[channel]['wavelength']), len(self.wavelength_response[channel]['response']))

        # Wavelength Response Interpolation could be done in this way
        # here more ions is good
        # response_interpolated =  interp1d(discrete_wavelengths,  self.wavelength_response[channel]['wavelength'], self.wavelength_response[channel]['response'])[0]
        # response_interpolated = np.interp(discrete_wavelengths,  self.wavelength_response[channel]['wavelength'], self.wavelength_response[channel]['response'])[0]

        # move to test: check: these are the same length because wavelength range in defining response functions
        # print('response', response, len(response), len(response[0]))
        # print(len(wrange))

        # print(len(temperature_array), len(temp_response), print(temp_response))
