import pytest

from sunpy.net.SPICE import attrs as a


def test_wrong_analysis_fk():
    with pytest.raises(ValueError,match = "given value must be boolean"):
        a.Analysis_fk(23)

def test_invalid_readme():
    value = 3
    with pytest.raises(ValueError,match = f"value must be boolean not {type(value)}"):
        a.Readme(value)

def test_invalid_kernel():
    with pytest.raises(ValueError,match="kernel type is required"):
        a.Kernel_type(None)
