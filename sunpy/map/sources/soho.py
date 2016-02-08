"""SOHO Map subclass definitions"""
from __future__ import absolute_import, print_function, division
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
from matplotlib import colors

from astropy.units import Quantity
from astropy.visualization import PowerStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.map import GenericMap
from sunpy.sun import constants
from sunpy.sun import sun
from sunpy.cm import cm

__all__ = ['EITMap', 'LASCOMap', 'MDIMap']


def _dsunAtSoho(date, rad_d, rad_1au=None):
    """Determines the distance to the Sun from SOhO following
    d_{\sun,Object} =
            D_{\sun\earth} \frac{\tan(radius_{1au}[rad])}{\tan(radius_{d}[rad])}
    though tan x ~ x for x << 1
    d_{\sun,Object} =
            D_{\sun\eart} \frac{radius_{1au}[rad]}{radius_{d}[rad]}
    since radius_{1au} and radius_{d} are dividing each other we can use [arcsec]
    instead.

    ---
    TODO: Does this apply just to observations on the same Earth-Sun line?
    If not it can be moved outside here.
    """
    if not rad_1au:
        rad_1au = sun.solar_semidiameter_angular_size(date)
    dsun = sun.sunearth_distance(date) * constants.au * (rad_1au / rad_d)
    # return scalar value not astropy.quantity
    return dsun.value


class EITMap(GenericMap):
    """SOHO EIT Image Map.

    SOHO EIT is an extreme ultraviolet (EUV) imager able to image the solar
    transition region and inner corona in four selected bandpasses,
    171 (Fe IX/X), 195 (Fe XII), 284 (Fe XV), and 304 (He II) Angstrom.

    SOHO was launched on 2 December 2 1995 into a sun-synchronous orbit and
    primary mission operations for SOHO EIT ended at the end of July 2010.

    References
    ----------
    * `SOHO Mission Page <http://sohowww.nascom.nasa.gov>`_
    * `SOHO EIT Instrument Page <http://umbra.nascom.nasa.gov/eit/>`_
    * `SOHO EIT User Guide <http://umbra.nascom.nasa.gov/eit/eit_guide/>`_
    """

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        # Fill in some missing info
        self.meta['detector'] = "EIT"
        self.meta['waveunit'] = "Angstrom"
        self._fix_dsun()
        self._nickname = self.detector
        self.plot_settings['cmap'] = cm.get_cmap(self._get_cmap_name())
        self.plot_settings['norm'] = ImageNormalize(stretch=PowerStretch(0.5))

    @property
    def rsun_obs(self):
        """
        Returns the solar radius as measured by EIT in arcseconds.
        """
        return Quantity(self.meta['solar_r'] * self.meta['cdelt1'], 'arcsec')

    def _fix_dsun(self):
        self.meta['dsun_obs'] = _dsunAtSoho(self.date, self.rsun_obs)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an EIT image"""
        return header.get('instrume') == 'EIT'


class LASCOMap(GenericMap):
    """SOHO LASCO Image Map

    The Large Angle and Spectrometric COronagraph (LASCO) is a set of three
    Lyot-type coronagraphs (C1, C2, and C3) that image the solar corona from
    1.1 to 32 solar radii.

    The C1 images rom 1.1 to 3 solar radii. The C2 telescope images the corona
    from 2 to 6 solar radii, overlaping the outer field-of-view of C1 from 2 to
    3 solar radii. The C3 telescope extends the field-of-view to 32 solar radii.

    SOHO was launched on 2 December 2 1995 into a sun-synchronous orbit.

    References
    ----------
    * `SOHO Mission Page <http://sohowww.nascom.nasa.gov>`_
    * `SOHO LASCO Instrument Page <http://lasco-www.nrl.navy.mil>`_
    * `SOHO LASCO Fits Header keywords <http://lasco-www.nrl.navy.mil/index.php?p=content/keywords>`_
    * `SOHO LASCO User Guide <http://lasco-www.nrl.navy.mil/index.php?p=content/handbook/hndbk>`_
    """

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        self.meta['CUNIT1'] = self.meta['CUNIT1'].lower()
        self.meta['CUNIT2'] = self.meta['CUNIT2'].lower()

        # Fill in some missing or broken info
        datestr = "{date}T{time}".format(date=self.meta.get('date-obs',
                                                            self.meta.get('date_obs')
                                                            ),
                                         time=self.meta.get('time-obs',
                                                            self.meta.get('time_obs')
                                                            )
                                         )
        self.meta['date-obs'] = datestr

        # If non-standard Keyword is present, correct it too, for compatibility.
        if 'date_obs' in self.meta:
            self.meta['date_obs'] = self.meta['date-obs']
        self.meta['wavelnth'] = np.nan
        self.meta['waveunit'] = 'nm'
        self._nickname = self.instrument + "-" + self.detector
        self.plot_settings['cmap'] = cm.get_cmap('soholasco{det!s}'.format(det=self.detector[1]))
        self.plot_settings['norm'] = ImageNormalize(stretch=PowerStretch(0.5))

    @property
    def measurement(self):
        """
        Returns the type of data taken.
        """
        # TODO: This needs to do more than white-light.  Should give B, pB, etc.
        return "white-light"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an LASCO image."""
        return header.get('instrume') == 'LASCO'


class MDIMap(GenericMap):
    """
    SOHO MDI Image Map

    The Michelson Doppler Imager (MDI) is a white light refracting telescope
    which feeds sunlight through a series of filters onto a CCD camera. Two
    tunable Michelson interformeters define a 94 mAngstrom bandpass that can be
    tuned across the Ni 6768 Angstrom solar absorption line.

    MDI measures line-of-sight motion (Dopplergrams), magnetic field
    (magnetograms), and brightness images of the full solar disk at several
    resolutions (4 arc-second to very low resolution) and a fixed selected
    region in higher resolution (1.2 arc-second).

    SOHO was launched on 2 December 2 1995 into a sun-synchronous orbit and
    SOHO MDI ceased normal science observations on 12 April 2011.

    References
    ----------
    * `SOHO Mission Page <http://sohowww.nascom.nasa.gov>`_
    * `SOHO MDI Instrument Page <http://soi.stanford.edu>`_
    * `SOHO MDI Fits Header keywords <http://soi.stanford.edu/sssc/doc/keywords.html>`_
    * `SOHO MDI Instrument Paper <http://soi.stanford.edu/sssc/doc/SP_paper_1995/MDI_SP_paper_1995.pdf>`_
    """

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        # Fill in some missing or broken info
        self.meta['detector'] = "MDI"
        self._fix_dsun()
        self.meta['wavelnth'] = np.nan
        self.meta['waveunit'] = 'nm'
        self._nickname = self.detector + " " + self.measurement
        vmin = np.nanmin(self.data)
        vmax = np.nanmax(self.data)
        if abs(vmin) > abs(vmax):
            self.plot_settings['norm'] = colors.Normalize(-vmin, vmin)
        else:
            self.plot_settings['norm'] = colors.Normalize(-vmax, vmax)


    @property
    def measurement(self):
        """
        Returns the type of data in the map.
        """
        return "magnetogram" if self.meta.get('content', " ").find('Mag') != -1 else "continuum"

    def _fix_dsun(self):
        """ Solar radius in arc-seconds at 1 au
            previous value radius_1au = 959.644
            radius = constants.average_angular_size
            There are differences in the keywords in the test FITS data and in
            the Helioviewer JPEG2000 files.  In both files, MDI stores the
            the radius of the Sun in image pixels, and a pixel scale size.
            The names of these keywords are different in the FITS versus the
            JP2 file.  The code below first looks for the keywords relevant to
            a FITS file, and then a JPEG2000 file.  For more information on
            MDI FITS header keywords please go to http://soi.stanford.edu/,
            http://soi.stanford.edu/data/ and
            http://soi.stanford.edu/magnetic/Lev1.8/ .
        """
        scale = self.meta.get('xscale', self.meta.get('cdelt1'))
        radius_in_pixels = self.meta.get('r_sun', self.meta.get('radius'))
        radius = scale * radius_in_pixels
        self.meta['radius'] = radius

        if not radius:
            # radius = sun.angular_size(self.date)
            self.meta['dsun_obs'] = constants.au
        else:
            self.meta['dsun_obs'] = _dsunAtSoho(self.date, radius)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an MDI image"""
        return header.get('instrume') == 'MDI' or header.get('camera') == 'MDI'
