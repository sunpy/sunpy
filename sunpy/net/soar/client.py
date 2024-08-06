"""
This file defines the SOARClient class which is used to access the Solar
Orbiter Archive (SOAR).
"""

import json
import pathlib
import re
from copy import copy
from json.decoder import JSONDecodeError

import astropy.table
import astropy.units as u
import requests
import sunpy.net.attrs as a
from sunpy import log
from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.time import parse_time

from sunpy.net.soar.attrs import SOOP, Product, walker

__all__ = ["SOARClient"]


class SOARClient(BaseClient):
    """
    Provides access to Solar Orbiter Archive (SOAR) which provides data for
    Solar Orbiter.

    References
    ----------
    * `SOAR <https://soar.esac.esa.int/soar/>`__
    """

    def search(self, *query, **kwargs):  # NOQA: ARG002
        r"""
        Query this client for a list of results.

        Parameters
        ----------
        *args: `tuple`
            `sunpy.net.attrs` objects representing the query.
        **kwargs: `dict`
            Any extra keywords to refine the search.
            Unused by this client.

        Returns
        -------
        A ``QueryResponseTable`` instance containing the query result.
        """
        query = and_(*query)
        queries = walker.create(query)

        results = []
        for query_parameters in queries:
            if "provider='SOAR'" in query_parameters:
                query_parameters.remove("provider='SOAR'")
            results.append(self._do_search(query_parameters))
        table = astropy.table.vstack(results)
        qrt = QueryResponseTable(table, client=self)
        qrt["Filesize"] = (qrt["Filesize"] * u.byte).to(u.Mbyte).round(3)
        qrt.hide_keys = ["Data item ID", "Filename"]
        return qrt

    @staticmethod
    def add_join_to_query(query: list[str], data_table: str, instrument_table: str):
        """
        Construct the WHERE, FROM, and SELECT parts of the ADQL query.

        Parameters
        ----------
        query : list[str]
            List of query items.
        data_table : str
            Name of the data table.
        instrument_table : str
            Name of the instrument table.

        Returns
        -------
        tuple[str, str, str]
            WHERE, FROM, and SELECT parts of the query.
        """
        final_query = ""
        # Extract wavemin and wavemax individually
        wavemin_pattern = re.compile(r"Wavemin='(\d+\.\d+)'")
        wavemax_pattern = re.compile(r"Wavemax='(\d+\.\d+)'")
        for current_parameter in query:
            parameter = copy(current_parameter)
            wavemin_match = wavemin_pattern.search(parameter)
            wavemax_match = wavemax_pattern.search(parameter)
            # If the wavemin and wavemax are same that means only one wavelength is given in query.
            if wavemin_match and wavemax_match and float(wavemin_match.group(1)) == float(wavemax_match.group(1)):
                # For PHI and SPICE, we can specify wavemin and wavemax in the query and get the results.
                # For PHI we have wavelength data in both angstrom and nanometer without it being mentioned in the SOAR.
                # For SPICE we get data in form of wavemin/wavemax columns, but only for the first spectral window.
                # To make sure this data is not misleading to the user we do not return any values for PHI AND SPICE.
                parameter = f"Wavelength='{wavemin_match.group(1)}'"
            elif wavemin_match and wavemax_match:
                parameter = f"Wavemin='{wavemin_match.group(1)}'+AND+h2.Wavemax='{wavemax_match.group(1)}'"
            prefix = "h1." if not parameter.startswith("Detector") and not parameter.startswith("Wave") else "h2."
            if parameter.startswith("begin_time"):
                time_list = parameter.split("+AND+")
                final_query += f"h1.{time_list[0]}+AND+h1.{time_list[1]}+AND+"
                # As there are no dimensions in STIX, the dimension index need not be included in the query for STIX.
                if "stx" not in instrument_table:
                    # To avoid duplicate rows in the output table, the dimension index is set to 1.
                    final_query += "h2.dimension_index='1'+AND+"
            else:
                final_query += f"{prefix}{parameter}+AND+"

        where_part = final_query[:-5]
        from_part = f"{data_table} AS h1"
        select_part = (
            "h1.instrument, h1.descriptor, h1.level, h1.begin_time, h1.end_time, "
            "h1.data_item_id, h1.filesize, h1.filename, h1.soop_name"
        )
        if instrument_table:
            from_part += f" JOIN {instrument_table} AS h2 USING (data_item_oid)"
            select_part += ", h2.detector, h2.wavelength, h2.dimension_index"
        return where_part, from_part, select_part

    @staticmethod
    def _construct_payload(query):
        """
        Construct search payload.

        Parameters
        ----------
        query : list[str]
            List of query items.

        Returns
        -------
        dict
            Payload dictionary to be sent with the query.
        """
        # Default data table
        data_table = "v_sc_data_item"
        instrument_table = None
        # Mapping is established between the SOAR instrument names and its corresponding SOAR instrument table alias.
        instrument_mapping = {
            "SOLOHI": "SHI",
            "EUI": "EUI",
            "STIX": "STX",
            "SPICE": "SPI",
            "PHI": "PHI",
            "METIS": "MET",
        }

        instrument_name = None
        for q in query:
            if q.startswith("instrument") or q.startswith("descriptor") and not instrument_name:
                instrument_name = q.split("=")[1][1:-1].split("-")[0].upper()
            elif q.startswith("level") and q.split("=")[1][1:3] == "LL":
                data_table = "v_ll_data_item"

        if instrument_name:
            if instrument_name in instrument_mapping:
                instrument_name = instrument_mapping[instrument_name]
            instrument_table = f"v_{instrument_name.lower()}_sc_fits"
            if data_table == "v_ll_data_item" and instrument_table:
                instrument_table = instrument_table.replace("_sc_", "_ll_")

        # Need to establish join for remote sensing instruments as they have instrument tables in SOAR.
        if instrument_name in ["EUI", "MET", "SPI", "PHI", "SHI"]:
            where_part, from_part, select_part = SOARClient.add_join_to_query(query, data_table, instrument_table)
        else:
            from_part = data_table
            select_part = "*"
            where_part = "+AND+".join(query)

        adql_query = {"SELECT": select_part, "FROM": from_part, "WHERE": where_part}

        adql_query_str = "+".join([f"{key}+{value}" for key, value in adql_query.items()])
        return {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "json", "QUERY": adql_query_str}

    @staticmethod
    def _do_search(query):
        """
        Query the SOAR server with a single query.

        Parameters
        ----------
        query : list[str]
            List of query items.

        Returns
        -------
        astropy.table.QTable
            Query results.
        """
        tap_endpoint = "http://soar.esac.esa.int/soar-sl-tap/tap"
        payload = SOARClient._construct_payload(query)
        # Need to force requests to not form-encode the parameters
        payload = "&".join([f"{key}={val}" for key, val in payload.items()])
        # Get request info
        r = requests.get(f"{tap_endpoint}/sync", params=payload, timeout=60)
        log.debug(f"Sent query: {r.url}")
        r.raise_for_status()
        try:
            response_json = r.json()
        except JSONDecodeError as err:
            msg = "The SOAR server returned an invalid JSON response. It may be down or not functioning correctly."
            raise RuntimeError(msg) from err

        # Do some list/dict wrangling
        names = [m["name"] for m in response_json["metadata"]]
        info = {name: [] for name in names}

        for entry in response_json["data"]:
            for i, name in enumerate(names):
                info[name].append(entry[i])

        if len(info["begin_time"]):
            info["begin_time"] = parse_time(info["begin_time"]).iso
            info["end_time"] = parse_time(info["end_time"]).iso

        result_table = astropy.table.QTable(
            {
                "Instrument": info["instrument"],
                "Data product": info["descriptor"],
                "Level": info["level"],
                "Start time": info["begin_time"],
                "End time": info["end_time"],
                "Data item ID": info["data_item_id"],
                "Filename": info["filename"],
                "Filesize": info["filesize"],
                "SOOP Name": info["soop_name"],
            },
        )
        if "detector" in info:
            result_table["Detector"] = info["detector"]
        if "wavelength" in info:
            result_table["Wavelength"] = info["wavelength"]
        result_table.sort("Start time")
        return result_table

    def fetch(self, query_results, *, path, downloader, **kwargs) -> None:  # NOQA: ARG002
        """
        Queue a set of results to be downloaded.
        `sunpy.net.base_client.BaseClient` does the actual downloading, so we
        just have to queue up the ``downloader``.

        Parameters
        ----------
        query_results : sunpy.net.fido_factory.UnifiedResponse
            Results from a Fido search.
        path : str
            Path to download files to. Must be a format string with a ``file``
            field for the filename.
        downloader : parfive.Downloader
            Downloader instance used to download data.
        kwargs :
            Keyword arguments aren't used by this client.
        """
        base_url = "http://soar.esac.esa.int/soar-sl-tap/data?" "retrieval_type=LAST_PRODUCT"

        for row in query_results:
            url = base_url
            if row["Level"].startswith("LL"):
                url += "&product_type=LOW_LATENCY"
            else:
                url += "&product_type=SCIENCE"
            data_id = row["Data item ID"]
            url += f"&data_item_id={data_id}"
            filepath = str(path).format(file=row["Filename"], **row.response_block_map)
            log.debug(f"Queuing URL: {url}")
            downloader.enqueue_file(url, filename=filepath)

    @classmethod
    def _can_handle_query(cls, *query) -> bool:
        """
        Check if this client can handle a given Fido query. Checks to see if a
        SOAR instrument or product is provided in the query.

        Returns
        -------
        bool
            True if this client can handle the given query.
        """
        required = {a.Time}
        optional = {a.Instrument, a.Detector, a.Wavelength, a.Level, a.Provider, Product, SOOP}
        if not cls.check_attr_types_in_query(query, required, optional):
            return False
        # check to make sure the instrument attr passed is one provided by the SOAR.
        # also check to make sure that the provider passed is the SOAR for which this client can handle.
        instr = [i[0].lower() for i in cls.register_values()[a.Instrument]]
        for x in query:
            if isinstance(x, a.Instrument) and str(x.value).lower() not in instr:
                return False
            if isinstance(x, a.Provider) and str(x.value).lower() != "soar":
                return False
        return True

    @classmethod
    def _attrs_module(cls):
        # Register SOAR specific attributes with Fido
        return "soar", "sunpy.net.soar.attrs"

    @classmethod
    def register_values(cls):
        """
        Register the SOAR specific attributes with Fido.

        Returns
        -------
        dict
            The dictionary containing the values formed into attributes.
        """
        return cls.load_dataset_values()

    @staticmethod
    def load_dataset_values():
        """
        Loads the net attribute values from the JSON file.

        Returns
        -------
        dict
            The dictionary containing the values formed into attributes.
        """
        # Instrument attrs
        attrs_path = pathlib.Path(__file__).parent / "data" / "attrs.json"
        with attrs_path.open() as attrs_file:
            all_datasets = json.load(attrs_file)
        # Convert from dict to list of tuples
        all_datasets = list(all_datasets.items())

        # Instrument attrs
        instr_path = pathlib.Path(__file__).parent / "data" / "instrument_attrs.json"
        with instr_path.open() as instr_attrs_file:
            all_instr = json.load(instr_attrs_file)
        all_instr = list(all_instr.items())

        soop_path = pathlib.Path(__file__).parent / "data" / "soop_attrs.json"
        with soop_path.open() as soop_path_file:
            all_soops = json.load(soop_path_file)

        all_soops = list(all_soops.items())

        return {
            Product: all_datasets,
            a.Instrument: all_instr,
            SOOP: all_soops,
            a.Provider: [("SOAR", "Solar Orbiter Archive.")],
        }
