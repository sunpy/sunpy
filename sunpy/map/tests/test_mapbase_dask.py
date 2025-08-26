import copy

import numpy as np
import pytest

import astropy.units as u
import astropy.wcs
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.wcs import FITSFixedWarning

import sunpy.map
from sunpy.data.test import get_test_filepath


@pytest.fixture
def aia171_test_dask_map(aia171_test_map):
    dask_array = pytest.importorskip('dask.array')
    return aia171_test_map._new_instance(
        dask_array.from_array(aia171_test_map.data),
        copy.deepcopy(aia171_test_map.meta)
    )


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_dask_array_repr(aia171_test_dask_map):
    # Check that _repr_html_ functions for a dask array
    html_dask_repr = aia171_test_dask_map._repr_html_(compute_dask=False)
    html_computed_repr = aia171_test_dask_map._repr_html_(compute_dask=True)
    assert html_dask_repr != html_computed_repr


# This is needed for the reproject_to function
with pytest.warns(VerifyWarning, match="Invalid 'BLANK' keyword in header."):  # NOQA: PT031
    with fits.open(get_test_filepath('aia_171_level1.fits')) as hdu:
        with pytest.warns(FITSFixedWarning, match="'datfix' made the change"):
            aia_wcs = astropy.wcs.WCS(header=hdu[0].header)


@pytest.mark.parametrize(("func", "args"), [
    ("max", {}),
    ("mean", {}),
    ("min", {}),
    pytest.param("reproject_to", {"wcs": aia_wcs}, marks=pytest.mark.xfail(reason="reproject is not dask aware")),
    pytest.param("resample", {"dimensions": (100, 100)*u.pix}, marks=pytest.mark.xfail()),
    pytest.param("rotate", {}, marks=pytest.mark.xfail(reason="nanmedian is not implemented in Dask")),
    ("std", {}),
    ("superpixel", {"dimensions": (10, 10)*u.pix}),
    ("submap", {"bottom_left": (100, 100)*u.pixel, "width": 10*u.pixel, "height": 10*u.pixel}),
])
def test_method_preserves_dask_array(aia171_test_map, aia171_test_dask_map, func, args):
    """
    Check that map methods preserve dask arrays if they are given as input, instead of eagerly
    operating on them and bringing them into memory.
    """
    # Check that result array is still a Dask array
    res_dask = aia171_test_dask_map.__getattribute__(func)(**args)
    res = aia171_test_map.__getattribute__(func)(**args)
    result_is_map = isinstance(res, sunpy.map.GenericMap)
    if result_is_map:
        assert isinstance(res_dask.data, type(aia171_test_dask_map.data))
        assert np.allclose(res_dask.data.compute(), res.data, atol=0.0, rtol=0.0)
    else:
        assert isinstance(res_dask, type(aia171_test_dask_map.data))
        assert np.allclose(res_dask.compute(), res, atol=0.0, rtol=0.0)
