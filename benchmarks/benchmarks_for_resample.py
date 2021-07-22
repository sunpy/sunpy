import astropy.units as u

import sunpy.data.sample
import sunpy.map


class TimeSuite:
    def setup(self):
        self._map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    def time_resample(self):
        self._map.resample([100, 100] * u.pix)


class PeakMemorySuite:
    def setup(self):
        self._map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    def peakmem_resample(self):
        self._map.resample([100, 100] * u.pix)
