import astropy.units as u
import sunpy.data.sample
import sunpy.map


class Creation:
    """
    Benchmarks for creating SunPy map objects.
    """
    params = ['AIA_171_IMAGE', 'HMI_LOS_IMAGE']
    param_names = ['name']

    def setup(self, name):
        """
        Set up the filename for the map creation based on the provided name.

        Args:
            name (str): The attribute name from `sunpy.data.sample` corresponding to a sample image.
        """
        self.filename = getattr(sunpy.data.sample, name)

    def time_create_map(self, name):
        """
        Benchmark the time taken to create a SunPy map.
        """
        sunpy.map.Map(self.filename)

    def mem_create_map(self, name):
        """
        Benchmark the memory usage when creating a SunPy map.

        Returns:
            Map: The created SunPy map.
        """
        return sunpy.map.Map(self.filename)

    def peakmem_create_map(self, name):
        """
        Benchmark the peak memory usage when creating a SunPy map.
        """
        sunpy.map.Map(self.filename)


class Resample:
    """
    Benchmarks for resampling SunPy maps.
    """

    def setup_cache(self):
        """
        Cache the SunPy map for resampling benchmarks.

        Returns:
            Map: The cached SunPy map.
        """
        return sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    def time_resample(self, aiamap):
        """
        Benchmark the time taken to resample a SunPy map.

        Args:
            aiamap (Map): The SunPy map to be resampled.
        """
        aiamap.resample([100, 100] * u.pix)

    def peakmem_resample(self, aiamap):
        """
        Benchmark the peak memory usage for resampling a SunPy map.

        Args:
            aiamap (Map): The SunPy map to be resampled.
        """
        aiamap.resample([100, 100] * u.pix)


class Rotate:
    """
    Benchmarks for rotating SunPy maps.
    """
    params = (['scipy', 'scikit-image', 'opencv'], range(0, 6))
    param_names = ['method', 'order']

    def setup_cache(self):
        """
        Cache the SunPy map for rotation benchmarks.

        Returns:
            Map: The cached SunPy map.
        """
        return sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    def setup(self, aiamap, method, order):
        """
        Validate the input parameters for the rotation benchmark.

        Args:
            aiamap (Map): The SunPy map to be rotated.
            method (str): The interpolation method.
            order (int): The order of the interpolation.

        Raises:
            NotImplementedError: If the method is 'opencv' and the order is unsupported.
        """
        if method == 'opencv' and order not in {0, 1, 3}:
            raise NotImplementedError("OpenCV only supports interpolation orders 0, 1, and 3.")

    def time_rotate(self, aiamap, method, order):
        """
        Benchmark the time taken to rotate a SunPy map.

        Args:
            aiamap (Map): The SunPy map to be rotated.
            method (str): The interpolation method.
            order (int): The order of the interpolation.
        """
        aiamap.rotate(30 * u.deg, method=method, order=order)

    def peakmem_rotate(self, aiamap, method, order):
        """
        Benchmark the peak memory usage for rotating a SunPy map.

        Args:
            aiamap (Map): The SunPy map to be rotated.
            method (str): The interpolation method.
            order (int): The order of the interpolation.
        """
        aiamap.rotate(30 * u.deg, method=method, order=order)
