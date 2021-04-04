"""
GOES SXI Map subclass definitions.
"""


from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.coordinates.sun import B0, angular_radius, earth_distance
from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch
from sunpy.sun import constants

__all__ = ['SXIMap']


class SXIMap(GenericMap):
    """
    GOES SXI Image Map

    The GOES Solar X-ray Imager (SXI) provides uninterrupted, full-disk, soft X-ray
    solar images, with a 1 min cadence and a single-image (adjustable) dynamic
    range near 100. A set of metallic thin-film filters provides temperature
    discrimination in the 0.6–6.0 nm bandpass. The spatial resolution of
    approximately 10 arcsec Full width half maximum (FWHM) is sampled with 5
    arcsec pixels.

    SXI was first launched on 23 July 2001 on NOAA’s GOES-12 satellite and has
    obtained uninterrupted images starting on 22 January 2003.

    .. note::

        The second image which is a pixel mask, is ignored and not read into the map.

    References
    ----------
    * `Mission Instrument Page <https://sxi.ngdc.noaa.gov/index.html>`__
    * `Data notes and Fits headers <https://sxi.ngdc.noaa.gov/sxi_data_notes.html>`__
    * `Official SXI Data Book <https://sxi.ngdc.noaa.gov/images/section06.pdf>`__
    * `Instrument Reference <https://doi.org/10.1007/s11207-005-7416-x>`__
    * `SXI Performance <https://doi.org/10.1007/s11207-005-7417-9>`__
    """
    # TODO: I added all of them, which might not be very useful.
    filter_wheel = ["OPEN", "PTHNA", "PTHNB", "PMEDB", "PMEDA", "PTHK",
                    "BTHNA", "BTHNB", "BMED", "BTHK", "UV", "RDSH"]

    def __init__(self, data, header, **kwargs):
        # Suggests that the header has been corrected already.
        if 'DETECTOR' not in header:
            # Check if the filter setting is in the list of expected values
            fw = header.get('WAVELNTH')
            if fw.lower().strip() not in [item.lower() for item in self.filter_wheel]:
                raise ValueError('Unexpected filter setting in header. '
                                 f'Expected {", ".join(self.filter_wheel)}. '
                                 f'Got {fw}')
            # Necessary fixes to the header
            header['DETECTOR'] = 'SXI'
            # Save the filter setting
            header['FILTER'] = header.get('WAVELNTH')
            # Update the wavelength with a number since a str is not supported.
            # TODO: change the number according to the filter
            header['WAVELNTH'] = "14"
            header['WAVEUNIT'] = 'Angstrom'
            # Need to switch these as they are not correct
            crota2 = header['CROTA1']
            crota1 = header['CROTA2']
            header['CROTA1'] = crota1
            header['CROTA2'] = crota2
            header['CTYPE1'] = "HPLN-TAN"
            header['CTYPE2'] = "HPLT-TAN"
            header['CUNIT1'] = "arcsec"
            header['CUNIT2'] = "arcsec"
            obs_time = header['date_obs']
            header['HGLT_OBS'] = B0(obs_time).to('deg').value
            header['HGLN_OBS'] = 0
            header['RSUN_OBS'] = angular_radius(obs_time).to('arcsec').value
            header['RSUN_REF'] = constants.radius.to('m').value
            header['DSUN_OBS'] = earth_distance(obs_time).to('m').value

        super().__init__(data, header, **kwargs)
        self._nickname = "{0}/{1}".format(self.meta['telescop'], self.detector)

        self.plot_settings['cmap'] = "sohoeit195"
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, AsinhStretch(0.01)), clip=False
        )

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """
        Determines if header corresponds to an GOES SXI image.
        """
        return str(header.get('instrume', '')).startswith('SXI')

    @property
    def measurement(self):
        """
        Returns the measurement type.
        """
        return self.meta['FILTER']
