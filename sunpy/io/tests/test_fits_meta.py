import pytest
import yaml

from astropy.table import Table

from sunpy.io.meta import fits_meta


def test_SolarHeader_class_basic_creation():
    """
    Test basic object creation
    """
    header = fits_meta.SolarnetHeader()
    assert isinstance(header, fits_meta.SolarnetHeader)


def test_loading_solarnet_schema():
    """
    Test load the solarnet schema file to find any errors in the file format
    """
    with open(fits_meta.solarnet_schema_file, "r") as f:
        fits_schema_yaml = yaml.safe_load(f)
    assert isinstance(fits_schema_yaml, dict)
    
    schema = fits_meta._schema_to_table(fits_schema_yaml)
    isinstance(schema, Table)

