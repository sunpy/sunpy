"""
Reads files in a .genx data directory.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import warnings
from collections import namedtuple

import h5py
import numpy as np
from astropy.table import QTable
import astropy.units as u
try:
    import ChiantiPy.core as ch
    import ChiantiPy.tools.data as ch_data
except ImportError:
    warnings.warn('Cannot import ChiantiPy. You will not be able to calculate the temperature response functions.')

from sunpy.io.special import read_genx

__author__ = ["Tessa D. Wilkinson", "Will Barnes"]


def aia_instr_properties_to_table(channel_list, instrument_files):
    """
    Read an AIA instrument (.genx) file into an Astropy table.

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


@u.quantity_input(temperature=u.K, density=u.cm**(-3), continuum_wavelength=u.angstrom)
def make_emiss_table(emiss_table_file, ion_list, temperature=np.logspace(5, 8, 50)*u.K,
                     density=1e15/np.logspace(5, 8, 50)*u.cm**(-3),
                     continuum_wavelength=np.arange(10., 500., 1.)*u.angstrom, **kwargs):
    """
    Build line and continuum emissivity table.

    Build HDF5 file of line and continuum contributions for each ion. For each ion, the
    contribution function as well as the three components of the continuum, two-photon,
    free-free, and free-bound, are calculated. Note that the contribution function is calculated
    here as,

    .. math::
       G(\lambda,n,T) = \\frac{1}{4\pi}\mathrm{Ab}(Z)\\frac{N(Z^{+z})}{N(Z)}\epsilon(\lambda,n,T)\\frac{1}{n_e},

    in units of photon :math:`\mathrm{cm}^{3}\,\mathrm{s}^{-1}\,\mathrm{str}^{-1}`. See the
    ChiantiPy documentation for more information regarding the continuum calculation.

    Parameters
    ----------
    emiss_table_file : `str`
    ion_list : `list`
        List of ion names, e.g. 'fe_15' for Fe XV
    temperature : array-like, optional
    density : array-like, optional
        Default is to assume a constant pressure of :math:`10^{15}\,\mathrm{K}\,\mathrm{cm}^{-3}`
    continuum_wavelength : array-like, optional

    See also
    --------
    EmissTableInterface : Easily access data in emissivity table files.
    """

    ch_data.Defaults['flux'] = 'photon'
    abundance_file = kwargs.get('abundance_file', 'sun_coronal_1992_feldman')
    ch_data.Defaults['ioneqfile'] = kwargs.get('ionization_equilibrium_file', ch_data.Defaults['ioneqfile'])

    with h5py.File(emiss_table_file, 'x') as hf:
        # save parameters
        dset_temperature = hf.create_dataset('temperature', data=temperature.value)
        dset_temperature.attrs['unit'] = temperature.unit.to_string()
        dset_density = hf.create_dataset('density', data=density.value)
        dset_density.attrs['unit'] = density.unit.to_string()
        dset_wavelength = hf.create_dataset('continuum_wavelength', data=continuum_wavelength.value)
        dset_wavelength.attrs['unit'] = continuum_wavelength.unit.to_string()

        for ion in ion_list:
            group = hf.create_group(ion)
            tmp_ion = ch.ion(ion, temperature=temperature.value, eDensity=density.value, abundance=abundance_file)
            # useful metadata
            group.attrs['spectroscopic_name'] = tmp_ion.Spectroscopic
            group.attrs['abundance_file'] = tmp_ion.AbundanceName
            group.attrs['ionization_equilibrium_file'] = tmp_ion.IoneqName
            group.attrs['Z'] = tmp_ion.Z
            group.attrs['stage'] = tmp_ion.Ion
            # line emission
            tmp_ion.emiss()
            emissivity = tmp_ion.Emiss['emiss'][np.argsort(tmp_ion.Emiss['wvl']), :].T*u.photon/u.s/u.steradian
            wavelength = np.sort(tmp_ion.Emiss['wvl'])*u.angstrom
            gofnt = tmp_ion.Abundance*emissivity*tmp_ion.IoneqOne[:,np.newaxis]/density[:,np.newaxis]
            # two-photon continuum
            tmp_ion.twoPhoton(continuum_wavelength.value)
            if 'rate' in tmp_ion.TwoPhoton:
                two_photon_continuum = tmp_ion.TwoPhoton['rate']
            else:
                two_photon_continuum = np.zeros((temperature.size, continuum_wavelength.size))
            # free-free continuum
            tmp_continuum = ch.Continuum(ion, temperature.value, abundance=abundance_file)
            tmp_continuum.calculate_free_free_emission(continuum_wavelength.value)
            free_free_continuum = tmp_continuum.free_free_emission
            # free-bound continuum
            try:
                tmp_continuum.calculate_free_bound_emission(continuum_wavelength.value)
                free_bound_continuum = tmp_continuum.free_bound_emission
            except ValueError:
                # free-bound information is not available for all ions
                free_bound_continuum = np.zeros((temperature.size, continuum_wavelength.size))
            # store results
            dset_gofnt = group.create_dataset('contribution_function', data=gofnt.value)
            dset_gofnt.attrs['unit'] = gofnt.unit.to_string()
            dset_transitions = group.create_dataset('transitions', data=wavelength.value)
            dset_transitions.attrs['unit'] = wavelength.unit.to_string()
            dset_tp = group.create_dataset('two_photon_continuum', data=two_photon_continuum)
            dset_tp.attrs['unit'] = 'photon cm^3 s^-1 sr^-1 angstrom^-1'
            dset_ff = group.create_dataset('free_free_continuum', data=free_free_continuum)
            dset_ff.attrs['unit'] = 'photon cm^3 s^-1 sr^-1 angstrom^-1'
            dset_fb = group.create_dataset('free_bound_continuum', data=free_bound_continuum)
            dset_fb.attrs['unit'] = 'photon cm^3 s^-1 sr^-1 angstrom^-1'


class EmissTableInterface(object):
    """
    Interface to the ion emission table.
    """

    def __init__(self,emiss_table_file):
        self.emiss_table_file = emiss_table_file

    @property
    def temperature(self):
        with h5py.File(self.emiss_table_file,'r') as hf:
            temperature = np.array(hf['temperature'])*u.Unit(hf['temperature'].attrs['unit'])
        return temperature

    @property
    def density(self):
        with h5py.File(self.emiss_table_file,'r') as hf:
            density = np.array(hf['density'])*u.Unit(hf['density'].attrs['unit'])
        return density

    @property
    def continuum_wavelength(self):
        with h5py.File(self.emiss_table_file,'r') as hf:
            continuum_wavelength = (np.array(hf['continuum_wavelength'])
                                    * u.Unit(hf['continuum_wavelength'].attrs['unit']))
        return continuum_wavelength

    def __getitem__(self, key):
        with h5py.File(self.emiss_table_file, 'r') as hf:
            if key not in hf:
                raise KeyError('{} not found in emission table {}'.format(key,self.emiss_table_file))
            ion_dict = dict(hf[key].attrs)
            for ds in hf[key]:
                ion_dict[ds] = np.array(hf['/'.join([key, ds])])*u.Unit(hf['/'.join([key, ds])].attrs['unit'])

        ion = namedtuple('ion', ' '.join([k for k in ion_dict.keys()]))
        return ion(**ion_dict)
