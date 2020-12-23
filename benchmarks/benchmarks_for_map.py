
import sunpy.data.sample
import sunpy.map


class TimeSuite:
    """
    An example benchmark that times the performance of sunpy map
    creation
    """
    def setup(self):
        self._map = sunpy.data.sample.AIA_171_IMAGE

    def time_create_sunpy_map(self):
        sunpy.map.Map(self._map)


class MemSuite:
    def setup(self):
        self._map = sunpy.data.sample.AIA_171_IMAGE

    def mem_sunpy_map(self):
        return sunpy.map.Map(self._map)


class PeakMemorySuite:
    def setup(self):
        self._map = sunpy.data.sample.AIA_171_IMAGE

    def peakmem_sunpy_map(self):
        sunpy.map.Map(self._map)
