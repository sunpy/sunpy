from sunpy.net import attrs
from sunpy.net.dataretriever import GenericClient
from sunpy.net.dataretriever.client import QueryResponse
from datetime import timedelta
import astropy.units as u


class AIASynopticClient(GenericClient):
    """
    A client for retrieving AIA synoptic data from JSOC.

    This client retrieves AIA synoptic data from the following source:
    https://jsoc1.stanford.edu/data/aia/synoptic/.
    """

    pattern = "https://jsoc1.stanford.edu/data/aia/synoptic/{{year:4d}}/{{year:4d}}{{month:2d}}{{day:2d}}_AIA_synoptic.fits"

    @classmethod
    def register_values(cls):
        return {
            attrs.Instrument: [("AIASynoptic", "AIA Synoptic data from the Solar Dynamics Observatory.")],
            attrs.Provider: [("JSOC", "Joint Science Operations Center.")],
        }

    @classmethod
    def _can_handle_query(cls, *query) -> bool:
        """Determines if the client can handle the given query."""
        required_attrs = {attrs.Time, attrs.Instrument}
        all_attrs = {type(attr) for attr in query}
        return required_attrs.issubset(all_attrs)

    def search(self, *query) -> QueryResponse:
        """Perform a search query for AIA synoptic data."""
        time_range = None

        for q in query:
            if isinstance(q, attrs.Time):
                time_range = q

        if not time_range:
            raise ValueError("Time range must be specified for the AIASynopticClient.")

        urls = self._generate_urls(time_range)
        return self._prepare_query_response(urls)

    def _generate_urls(self, time_range: attrs.Time) -> list:
        """Generate a list of URLs for AIA synoptic data."""
        urls = []
        current_time = time_range.start.datetime

        while current_time <= time_range.end.datetime:
            urls.append(
                self.pattern.format(year=current_time.year, month=current_time.month, day=current_time.day)
            )
            current_time += timedelta(days=1)  # Increment by a day

        return urls

    def _prepare_query_response(self, urls: list) -> QueryResponse:
        """Prepare the query response."""
        data = {
            "Start Time": [url.split("/")[-1].split("_")[0] for url in urls],
            "Instrument": ["AIASynoptic"] * len(urls),
            "url": urls,
        }
        return QueryResponse(data)


def main():
    start_time = "2020-01-01T00:00:00.00"
    end_time = "2020-01-01T23:59:59.999"

    # Client = AIASynopticClient()
    from sunpy.net import Fido, attrs as a

    # from sunpy.net import Fido
    # from sunpy.net.dataretriever.sources.aia_synoptic import AIASynopticClient

    # Register the AIASynopticClient with Fido
    # Fido.register_client(AIASynopticClient)

    request = Fido.search(
        attrs.Instrument("AIA"),
        attrs.Time(start_time, end_time),
        attrs.Sample(1 * u.hour),
        attrs.Wavelength(171 * u.angstrom),
    )
    print(request)


if __name__ == "__main__":
    main()
