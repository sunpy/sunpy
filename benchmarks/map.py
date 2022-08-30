import astropy.units as u

import sunpy.data.sample
import sunpy.map


class Creation:
    params = ['AIA_171_IMAGE', 'HMI_LOS_IMAGE']
    param_names = ['name']

    def setup(self, name):
        self.filename = getattr(sunpy.data.sample, name)

    def time_create_map(self):
        sunpy.map.Map(self.filename)

    def mem_create_map(self):
        return sunpy.map.Map(self.filename)

    def peakmem_create_map(self):
        sunpy.map.Map(self.filename)


class Resample:
    def setup_cache(self):
        aiamap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
        return aiamap

    def time_resample(self, aiamap):
        aiamap.resample([100, 100] * u.pix)

    def peakmem_resample(self, aiamap):
        aiamap.resample([100, 100] * u.pix)


class Rotate:
    params = (['scipy', 'scikit-image', 'opencv'], range(0, 6))
    param_names = ['method', 'order']

    def setup_cache(self):
        aiamap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
        return aiamap

    def time_rotate(self, aiamap, method, order):
        aiamap.rotate(30*u.deg, method=method, order=order)

    def peakmem_rotate(self, aiamap, method, order):
        aiamap.rotate(30*u.deg, method=method, order=order)
