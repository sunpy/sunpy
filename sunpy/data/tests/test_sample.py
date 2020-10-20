import pytest

from .._sample import _download_sample_data, _retry_sample_data


@pytest.mark.remote_data
def test_retry_sample_data(tmpdir):
    # Wrong base URL.
    result = _download_sample_data(
        "http://ipv4.download.thinkbroadband.com", [("tca110607.fits",
                                                     tmpdir.strpath+"/tca110607.fits")], False)
    assert result == []
    assert result.errors != []

    # Inside _retry_sample_data, the base url will be updated to point towards
    # http://data.sunpy.org/sunpy/v1/
    result_retry = _retry_sample_data(result)
    assert result_retry == [tmpdir.strpath+'/tca110607.fits']
    assert result_retry.errors == []


@pytest.mark.remote_data
def test_download_sample_data(tmpdir):
    # Download a simple random file off the internet.
    result = _download_sample_data(
        "http://ipv4.download.thinkbroadband.com", [("5MB.zip", tmpdir.strpath+"/5MB.zip")], False)
    assert result == [tmpdir.strpath+"/5MB.zip"]
    assert result.errors == []
