import os

import pytest
from hypothesis import given

import astropy.units as u

import sunpy.net.dataretriever.sources.aia_synoptic as aia_synoptic
from sunpy.net import attrs as a
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.tests.strategies import time_attr
from sunpy.time import TimeRange, parse_time


@pytest.fixture
def aia_synoptic_client():
    """Fixture to provide an instance of AIASynopticClient."""
    return aia_synoptic.AIASynopticClient()


@given(time_attr())
def test_can_handle_query(aia_synoptic_client, time):
    """
    Test if the AIASynopticClient can handle specific queries.
    """

    # Test with synoptic data
    assert aia_synoptic_client._can_handle_query(time, a.Instrument.aia, aia_synoptic.AIASynopticData())

    # If AIASynopticData is present, ignore resolution mismatch
    assert aia_synoptic_client._can_handle_query(
        time, a.Instrument.aia, a.Resolution(4.0), aia_synoptic.AIASynopticData()
    )

    # Test normal query without synoptic data
    assert aia_synoptic_client._can_handle_query(time, a.Instrument.aia)
    assert aia_synoptic_client._can_handle_query(time, a.Instrument.aia, a.Resolution(4.0))
    assert not aia_synoptic_client._can_handle_query(time)


def mock_query_object(
    aia_synoptic_client, start_time="2020-01-01T00:00:00.00", end_time="2020-01-01T23:59:59.999"
):
    """
    Creating a Query Response object and prefilling it with some information.

    Parameters:
        aia_synoptic_client (AIASynopticClient): Client instance used for query response.
        start_time (str): Start time for the mocked query.
        end_time (str): End time for the mocked query.

    Returns:
        QueryResponse: A mocked query response object.
    """
    obj = {
        "Start Time": parse_time(start_time),
        "End Time": parse_time(end_time),
        "Instrument": "AIASynoptic",
        "Physobs": "intensity",
        "Source": "SDO",
        "Provider": "JSOC",
        "Wavelength": "0171",
        "url": (r"https://jsoc1.stanford.edu/data/aia/synoptic/2023/10/11/H0000/AIA20231011_0000_0193.fits"),
    }
    return QueryResponse([obj], client=aia_synoptic_client)


@pytest.mark.remote_data
def test_fetch_working(aia_synoptic_client, tmp_path):
    """
    Tests that data is fetched correctly.

    Parameters:
        aia_synoptic_client (AIASynopticClient): Fixture client instance.
        tmp_path: Pytest fixture for creating temporary directories.
    """
    response = mock_query_object(aia_synoptic_client)
    download_dir = tmp_path / "downloads"

    # Ensure the directory is empty before downloading
    assert len(list(download_dir.glob("*"))) == 0

    files = aia_synoptic_client.fetch(response, path=str(download_dir))

    # Ensure exactly one file was downloaded
    assert len(files) == 1
    assert files[0].endswith(".fits")

    # Ensure the file is correctly downloaded
    downloaded_file = download_dir / os.path.basename(files[0])
    assert downloaded_file.exists()
    assert downloaded_file.stat().st_size > 0  # Check if file size is greater than 0


def test_basic_query(aia_synoptic_client):
    """
    Basic test to perform a search query for synoptic data.
    """
    # Define a time range for the query
    time_range = TimeRange("2020-01-01", "2020-01-02")

    # Perform a search query
    query_result = aia_synoptic_client.search(
        a.Time(time_range.start, time_range.end),
        a.Wavelength(171 * u.angstrom),
        a.Sample(12 * u.hour),
        a.Instrument("AIA"),
        aia_synoptic.AIASynopticData(),
    )

    # Output the results
    assert query_result is not None
    print(f"Query Result: {query_result}")

    # Fetch the data for the first result (optional)
    if query_result:
        download_path = "./downloads"
        files = aia_synoptic_client.fetch(query_result[:1], path=download_path, downloader=None)
        assert len(files) > 0
        print(f"Downloaded files: {files}")


if __name__ == "__main__":
    from sunpy.net.dataretriever.sources.aia_synoptic import AIASynopticClient
    from sunpy.net import attrs

    start_time = "2020-01-01T00:00:00.00"
    end_time = "2020-01-01T23:59:59.999"

    Client = AIASynopticClient()
    search_request = Client.search(
        attrs.Time(start_time, end_time), attrs.Sample(1 * u.hour), attrs.Wavelength(193 * u.angstrom)
    )
    request = Client.fetch(search_request)
    print(search_request)
    print(request)
