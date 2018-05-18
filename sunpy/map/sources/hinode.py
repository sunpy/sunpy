"""Hinode XRT and SOT Map subclass definitions"""
from __future__ import absolute_import, print_function, division
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = ["Jack Ireland, Jose Ivan Campos-Rozo, David Perez-Suarez"]
__email__ = "jack.ireland@nasa.gov"

import matplotlib.pyplot as plt

from sunpy.map import GenericMap

__all__ = ['XRTMap', 'SOTMap']

# the following values comes from xrt_prep.pro
# search for saturation in
# http://darts.jaxa.jp/pub/ssw/hinode/xrt/idl/util/xrt_prep.pro
# SATURATION_LIMIT = 2500


def _lower_list(l):
    return [item.lower() for item in l]


class XRTMap(GenericMap):
    """Hinode XRT map definition.

    The X-Ray Telescope (XRT) is a high resolution grazing incidence telescope,
    which is a succsessor to Yohkoh. It provides 2-arcsecond resolution images
    of the highest temperature solar coronal material,
    from 1,000,000 to 10,000,000 Kelvin.

    Hinode was launched on 22 September 2006 into a sun-synchronous orbit.

    References
    ----------
    * `Hinode Mission Page <https://solarb.msfc.nasa.gov/index.html>`_
    * `XRT Instrument Page <http://xrt.cfa.harvard.edu>`_
    * `Fits header reference <http://hinode.nao.ac.jp/uploads/2016/04/22/SB_MW_Key13.pdf>`_
    * `Hinode User Guide <http://hinode.nao.ac.jp/en/for-researchers/analysis-guide/>`_
    * `XRT Analysis Guide <http://xrt.cfa.harvard.edu/science/tutorials.php>`_
    * `Coronal Temperature Diagnostic Capability of the Hinode/X-Ray Telescope Based on Self-Consistent Calibration <https://arxiv.org/abs/1011.2867>`_
    """
    filter_wheel1_measurements = ["Al_med", "Al_poly", "Be_med",
                                  "Be_thin", "C_poly", "Open"]
    filter_wheel2_measurements = ["Open", "Al_mesh", "Al_thick",
                                  "Be_thick", "Gband", "Ti_poly"]

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        # converting data array to masked array
        # self.data = ma.masked_where(self.data > SATURATION_LIMIT, self.data)

        fw1 = header.get('EC_FW1_')
        if fw1.lower() not in _lower_list(self.filter_wheel1_measurements):
            raise ValueError('Unpexpected filter wheel 1 in header.')
        fw1 = fw1.replace("_", " ")

        fw2 = header.get('EC_FW2_')
        if fw2.lower() not in _lower_list(self.filter_wheel2_measurements):
            raise ValueError('Unpexpected filter wheel 2 in header.')
        fw2 = fw2.replace("_", " ")

        self.meta['detector'] = "XRT"
#        self.meta['instrume'] = "XRT"
        self.meta['telescop'] = "Hinode"
        self.plot_settings['cmap'] = plt.get_cmap(name='hinodexrt')

    @property
    def measurement(self):
        fw1 = self.meta.get('EC_FW1_').replace("_", " ")
        fw2 = self.meta.get('EC_FW2_').replace("_", " ")
        return "{0}-{1}".format(fw1, fw2)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an XRT image"""
        return header.get('instrume') == 'XRT'


class SOTMap(GenericMap):
    """Hinode SOT Image Map definition.

    The Hinode Solar Optical Telescope (SOT) consists of a 50 cm
    diffraction-limited Gregorian telescope. It is optimized for
    accurate measurement of the vector magnetic
    field in the photosphere and dynamics of both the photosphere and
    chromosphere associated with the magnetic fields.

    Hinode was launched on 22 September 2006 into a sun-synchronous orbit.

    References
    ----------
    * `Hinode Mission Page <http://solarb.msfc.nasa.gov/index.html>`_
    * `Hinode SOT Instrument Page <http://sot.lmsal.com>`_
    * `Hinode SOT Instrument Paper <https://arxiv.org/abs/0711.1715>`_
    * `Data Analsis Guide <https://sot.lmsal.com/doc/rep/sot254/fid366/SOT00042_33_SOT_Analysis_Guide_SAG.pdf>`_
    """
    # TODO: get a link for the SOT FITS headers
    # Add in some information about the the possible instrument, observation
    # type, observable ion and wavelength

    Instruments = ['SOT/WB', 'SOT/NB', 'SOT/SP', 'SOT/CT']

    Waves = ['6302A', 'BFI no move', 'CN bandhead 3883',
             'Ca II H line', 'G band 4305', 'NFI no move',
             'TF Fe I 6302', 'TF Mg I 5172', 'TF Na I 5896',
             'blue cont 4504', 'green cont 5550', 'red cont 6684']

    Observation_Type = ['FG (simple)', 'FG focus scan',
                        'FG shuttered I and V', 'FG shutterless I and V',
                        'FG shutterless I and V with 0.2s intervals',
                        'FG shutterless Stokes', 'SP IQUV 4D array']

    def __init__(self, data, header, **kwargs):
        GenericMap.__init__(self, data, header, **kwargs)

        self.meta['detector'] = "SOT"
        self.meta['telescop'] = "Hinode"
        self._nickname = self.detector

        # TODO (add other options, Now all threated as intensity. This follows
        # Hinode SDC archive) StokesQUV -> grey, Velocity -> EIS, Width -> EIS,
        # Mag Field Azi -> IDL 5 (STD gamma II)
        # 'WB' -> red
        # 'NB'(0 = red); (>0 = gray), # nb has 1 stokes I, the rest quv
        # 'SP' (<=1 = red); (>1 = gray) #sp has 2 stokes I, the rest quv
        color = {'SOT/WB': 'intensity',
                 'SOT/NB': 'intensity',  # For the 1st dimension
                 'SOT/SP': 'intensity',  # For the 1st 2 dimensions
                 }

        self.plot_settings['cmap'] = plt.get_cmap('hinodesot' + color[self.instrument])

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an SOT image."""
        return header.get('instrume') in cls.Instruments
