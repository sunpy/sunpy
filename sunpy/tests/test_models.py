import pytest
import time
import pathlib
from astropy.table import QTable
import sunpy.sun.models  as sun_models  

_MODEL_DATA_DIR = pathlib.Path(sun_models.__file__).parent / "data"


_MODELS = {
    "chromosphere_avrett_loeser_2008": _MODEL_DATA_DIR / "chromosphere_avrett_Loeser_2008_model.ecsv",
  
}

def test_valid_model_loading():
    """
    Test if models in _MODELS load correctly and are cached.
    """
    for model_name in _MODELS.keys():
        model_data = getattr(sun_models, model_name)  # Ensure attribute exists
        assert model_data is not None, f"Test Failed: {model_name} does not exist in sun_models."
        assert isinstance(model_data, QTable), f"Test Failed: {model_name} should return a QTable object."
        assert model_name in sun_models._MODEL_CACHE, f"Test Failed: {model_name} should be cached after loading."

def test_cache_functionality():
    """
    Ensure model caching works properly.
    """
    for model_name in _MODELS.keys():
        model_data_1 = getattr(sun_models, model_name)
        model_data_2 = getattr(sun_models, model_name)
        assert model_data_1 is model_data_2, f"Test Failed: Cache is not returning the same object for {model_name}."

def test_invalid_model_handling():
    """
    Ensure an invalid model raises AttributeError.
    """
    with pytest.raises(AttributeError):  
        getattr(sun_models, "invalid_model")  

def test_performance():
    """
    Ensure model loads within an acceptable time.
    """
    for model_name in _MODELS.keys():
        start_time = time.time()
        getattr(sun_models, model_name)  
        duration = time.time() - start_time
        assert duration < 1.0, f"Test Failed: Model '{model_name}' loading took too long! ({duration:.4f} seconds)"

# Run tests
if __name__ == "__main__":
    pytest.main(["-v"])  # Verbose output
