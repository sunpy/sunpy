from io import BytesIO

from asv_runner.benchmarks.mark import SkipNotImplemented, skip_benchmark
from matplotlib.figure import Figure

import astropy.units as u

import sunpy.data.sample
import sunpy.map
from sunpy.coordinates import propagate_with_solar_surface


class Creation:
    params = ['AIA_171_IMAGE', 'HMI_LOS_IMAGE']
    param_names = ['name']

    def setup(self, name):
        self.filename = getattr(sunpy.data.sample, name)

    def time_create_map(self, name):
        sunpy.map.Map(self.filename)

    # Skipped due to a bug in pympler.asizeof
    # https://github.com/pympler/pympler/issues/151
    @skip_benchmark
    def mem_create_map(self, name):
        return sunpy.map.Map(self.filename)

    def peakmem_create_map(self, name):
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

    def setup(self, aiamap, method, order):
        if method == 'opencv' and order not in {0, 1, 3}:
            raise SkipNotImplemented

    def time_rotate(self, aiamap, method, order):
        aiamap.rotate(30*u.deg, method=method, order=order)

    def peakmem_rotate(self, aiamap, method, order):
        aiamap.rotate(30*u.deg, method=method, order=order)


class Reproject:
    params = ['interpolation', 'adaptive']
    param_names = ['algorithm']

    def setup_cache(self):
        maps = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE, sunpy.data.sample.HMI_LOS_IMAGE)
        return maps

    def time_reproject_to(self, maps, algorithm):
        maps[1].reproject_to(maps[0].wcs, algorithm=algorithm)

    def peakmem_reproject_to(self, maps, algorithm):
        maps[1].reproject_to(maps[0].wcs, algorithm=algorithm)

    def time_reproject_to_plus_diffrot(self, maps, algorithm):
        with propagate_with_solar_surface():
            maps[1].reproject_to(maps[0].wcs, algorithm=algorithm)

    def peakmem_reproject_to_plus_diffrot(self, maps, algorithm):
        with propagate_with_solar_surface():
            maps[1].reproject_to(maps[0].wcs, algorithm=algorithm)


class Autoalign:
    params = [False, True, 'pcolormesh']
    param_names = ['autoalign']

    def setup_cache(self):
        maps = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE, sunpy.data.sample.HMI_LOS_IMAGE)
        return maps

    def time_autoalign(self, maps, autoalign):
        fig = Figure()
        ax = fig.add_subplot(projection=maps[0])
        maps[1].plot(axes=ax, autoalign=autoalign)
        buf = BytesIO()
        fig.savefig(buf, format='png')

    def peakmem_autoalign(self, maps, autoalign):
        fig = Figure()
        ax = fig.add_subplot(projection=maps[0])
        maps[1].plot(axes=ax, autoalign=autoalign)
        buf = BytesIO()
        fig.savefig(buf, format='png')
