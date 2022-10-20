import pytest

from sunpy.data._sample import _download_sample_data, _retry_sample_data


@pytest.mark.remote_data
def test_retry_sample_data(tmpdir):
    # Wrong base URL.
    result = _download_sample_data(
        "http://ipv4.download.thinkbroadband.com", [("tca110607.fits",
                                                     tmpdir / "tca110607.fits")], False)
    assert result == []
    assert result.errors != []

    result_retry = _retry_sample_data(result, "http://data.sunpy.org/sunpy/v1/")
    assert result_retry == [tmpdir / 'tca110607.fits']
    assert result_retry.errors == []


@pytest.mark.remote_data
def test_download_sample_data(tmpdir):
    # Download a simple random file off the internet.
    result = _download_sample_data(
        "http://ipv4.download.thinkbroadband.com", [("5MB.zip", tmpdir / "5MB.zip")], False)
    assert result == [tmpdir / "5MB.zip"]
    assert result.errors == []
