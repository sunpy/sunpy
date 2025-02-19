import pytest

import sunpy.sun.models as sun_models


@pytest.mark.parametrize("model_name", sun_models._MODELS.keys())
def test_valid_model_loading(model_name):
    """
    Test if models in sun_models._MODELS load correctly.
    """
    model_data = getattr(sun_models, model_name)
    assert model_data is not None


def test_cache_functionality():
    """
    Ensure model caching works properly for a specific model.
    """
    model_name = "avrett_loeser_2008"
    assert model_name in sun_models._MODELS

    model_data_1 = getattr(sun_models, model_name)
    model_data_2 = getattr(sun_models, model_name)

    assert model_data_1 is model_data_2
    assert model_name in sun_models._MODEL_CACHE


def test_invalid_model_handling():
    """
    Ensure an invalid model raises AttributeError with the correct message.
    """
    with pytest.raises(AttributeError, match="Error: Model 'invalid_model' is not available."):
        getattr(sun_models, "invalid_model")
