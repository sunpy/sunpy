import pytest
import sunpy.sun.models  as sun_models  

@pytest.mark.parametrize("model_name", sun_models._MODELS.keys())
def test_valid_model_loading(model_name):
    """
    Test if models in sun_models._MODELS load correctly.
    """
    model_data = getattr(sun_models, model_name)
    assert model_data is not None, f"Test Failed: {model_name} does not exist."
    

@pytest.mark.parametrize("model_name", sun_models._MODELS.keys())
def test_cache_functionality(model_name):
    """
    Ensure model caching works properly.
    """
    model_data_1 = getattr(sun_models, model_name)
    model_data_2 = getattr(sun_models, model_name)
    assert model_data_1 is model_data_2, f"Test Failed: Cache is not returning the same object for {model_name}."


def test_invalid_model_handling():
    """
    Ensure an invalid model raises AttributeError with the correct message.
    """
    with pytest.raises(AttributeError, match="Error: Model 'invalid_model' is not available."):
        getattr(sun_models, "invalid_model")  
 

if __name__ == "__main__":
    pytest.main(["-v"])  
