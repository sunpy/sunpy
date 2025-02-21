import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import Longitude
from astropy.tests.helper import assert_quantity_allclose

from sunpy.sun.models import differential_rotation
import sunpy.sun.models as sun_models


@pytest.fixture
def seconds_per_day():
    return 24 * 60 * 60.0 * u.s


def test_single(seconds_per_day):
    rot = differential_rotation(10 * seconds_per_day, 30 * u.deg)
    assert_quantity_allclose(rot, 136.8216 * u.deg, rtol=1e-3)


def test_array(seconds_per_day):
    rot = differential_rotation(10 * seconds_per_day, np.linspace(-70, 70, 2) * u.deg)
    assert_quantity_allclose(rot, Longitude(np.array([110.2725, 110.2725]) * u.deg), rtol=1e-3)


def test_synodic(seconds_per_day):
    rot = differential_rotation(10 * seconds_per_day, 30 * u.deg, model='howard', frame_time='synodic')
    assert_quantity_allclose(rot, 126.9656 * u.deg, rtol=1e-3)


def test_sidereal(seconds_per_day):
    rot = differential_rotation(10 * seconds_per_day, 30 * u.deg, model='howard', frame_time='sidereal')
    assert_quantity_allclose(rot, 136.8216 * u.deg, rtol=1e-3)


def test_howard(seconds_per_day):
    rot = differential_rotation(10 * seconds_per_day, 30 * u.deg, model='howard')
    assert_quantity_allclose(rot, 136.8216 * u.deg, rtol=1e-3)


def test_allen(seconds_per_day):
    rot = differential_rotation(10 * seconds_per_day, 30 * u.deg, model='allen')
    assert_quantity_allclose(rot, 136.9 * u.deg, rtol=1e-3)


def test_snodgrass(seconds_per_day):
    rot = differential_rotation(10 * seconds_per_day, 30 * u.deg, model='snodgrass')
    assert_quantity_allclose(rot, 135.4232 * u.deg, rtol=1e-3)


def test_rigid(seconds_per_day):
    rot = differential_rotation(10 * seconds_per_day, [0, 30, 60] * u.deg, model='rigid')
    assert_quantity_allclose(rot, [141.844 * u.deg] * 3, rtol=1e-3)


def test_fail(seconds_per_day):
    with pytest.raises(ValueError, match="model must equal one of { howard , snodgrass , allen , rigid }"):
        differential_rotation(10 * seconds_per_day, 30 * u.deg, model='garbage')


@pytest.mark.parametrize("model_name", sun_models._MODELS.keys())
def test_valid_model_loading(model_name):
    model_data = getattr(sun_models, model_name)
    assert model_data is not None


def test_cache_functionality():
    model_name = "chromosphere_avrett_loeser_2008"
    model_data_1 = getattr(sun_models, model_name)
    model_data_2 = getattr(sun_models, model_name)

    assert model_data_1 is model_data_2
    assert model_name in sun_models._MODEL_CACHE


def test_invalid_model_handling():
    with pytest.raises(AttributeError, match="Error: Model 'invalid_model' is not available."):
        getattr(sun_models, "invalid_model")
