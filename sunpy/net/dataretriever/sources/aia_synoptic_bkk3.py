from sunpy.net import attrs as a
from sunpy.net.dataretriever import GenericClient
import astropy.units as u

__all__ = ["AIASynopticClient"]


class AIASynopticClient(GenericClient):
    """
    A client for retrieving AIA synoptic data from JSOC.
    """

    baseurl = r"https://jsoc1.stanford.edu/data/aia/synoptic/%Y/%m/%d/H%H00/"
    pattern = "{}AIA{year}{month:02d}{day:02d}_{hour:02d}{minute:02d}_{Wavelength}{}"
    known_wavelengths = ["0094", "0131", "0171", "0193", "0211", "0304", "1600", "1700"]

    @classmethod
    def _can_handle_query(cls, *query) -> bool:
        required_attrs = {a.Time, a.Instrument}
        return required_attrs.issubset({type(attr) for attr in query})

    @classmethod
    def register_values(cls):
        adict = {
            a.Instrument: [("AIA", "Atmospheric Imaging Assembly"), ("aiasynoptic", "AIA Synoptic Series")],
            a.Wavelength: [[str(wl).zfill(4), f"{wl} Angstroms"] for wl in cls.known_wavelengths],
            # a.ExtentType: [("SYNOPTIC", "Coverage of a complete solar rotation synthesized over time")],
        }
        return adict

    @staticmethod
    def check_file_exists(url):
        response = requests.head(url)
        return response.status_code == 200


def main():
    from sunpy.net import Fido
    from sunpy.net import attrs as a
    import astropy.units as u

    # Define your query parameters
    time_range = a.Time("2023-10-11", "2023-10-12")
    instrument = a.Instrument("AIASynoptic")
    wavelength = a.Wavelength(94 * u.angstrom)

    # Perform the search
    request = Fido.search(time_range, instrument, wavelength)

    print("Found the following requests:")
    print(request)

    # Check if any results were returned
    if request:
        for resp in request.responses:
            if hasattr(resp, "table"):
                print("Column names:", resp.table.colnames)
            else:
                print("No table available in this response.")
    else:
        print("No results found for the query.")


if __name__ == "__main__":
    main()
