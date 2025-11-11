import os
import json
import urllib.request
from pathlib import Path

import requests

from sunpy.net import attrs as a
from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.net.solarnet.attrs import Dataset, Tags, Target, walker

_BASE_URL = "https://solarnet2.oma.be/service/api/svo/{}"

__all__ = ["SOLARNETClient"]

class SOLARNETClient(BaseClient):
    """
    Provides access to query and download from the SOLARNET Virtual Observatory (SVO).

    The SVO was first supported by the `SOLARNET <http://solarnet-east.eu/>`__ project,
    funded by the European Commission's FP7 Capacities Programme under the Grant Agreement 312495.
    It was turned operational due to the `SOLARNET2 <https://solarnet-project.eu/>`__ project, funded by the European Union's Horizon 2020
    Research and Innovation Programme under Grant Agreement 824135.
    It's purpose is to collect metadata from as many solar observations as possible,
    especially those involved in the SOLARNET projects, in a common catalog and make them available to the scientific community.

    There is various ways one can access SVO `Data <https://solarnet.oma.be/#introduction>`__, this client using the
    RESTful API (`available documentation <https://solarnet.oma.be/svo_restful_api_user_manual.html>`__).

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> search_results = Fido.search(a.solarnet.Dataset.gris_level_1,a.solarnet.Limit(2), a.solarnet.Target.ar)   # doctest: +REMOTE_DATA
    >>> search_results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the SOLARNETClient:
    Source: https://solarnet2.oma.be
    <BLANKLINE>
         Oid               Date_beg                 Date_end          Wavemin   Wavemax  Accumula Adunorm Aosystem Azimut B0angle Bscale Bzero Camera Campos Coefi_1 Coefi_2 Coefi_3 Coefi_4 Coefq_1 Coefq_2 Coefq_3 Coefq_4 Coefu_1 Coefu_2 Coefu_3 Coefu_4 Coefv_1 Coefv_2 Coefv_3 Coefv_4 Colpos Datafram Datavers         Date_obs         Dec  Eimgacc Elevatio Elperadu Exptime Filename Filestat Fullfram Gratangl Hourangl Ifile Imagtype Irfilter Iserie L0angle Lc1_1 Lc1_2 Lc1_3 Lc1_4 Lc2_1 Lc2_2 Lc2_3 Lc2_4 Lcs  Measure Naxis Naxis1 Naxis2 Naxis3 Nfiles P0angle Paraangl R0radius  Ra  Rotangle Rotcode Series Spsystem States Stepangl Steps Stepsh Stepsize Stepsv Telescop Temp_lc1 Temp_lc2 Waveleng F_size
    -------------- ------------------------ ------------------------ --------- --------- -------- ------- -------- ------ ------- ------ ----- ------ ------ ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------ -------- -------- ------------------------ ---- ------- -------- -------- ------- -------- -------- -------- -------- -------- ----- -------- -------- ------ ------- ----- ----- ----- ----- ----- ----- ----- ----- ---- ------- ----- ------ ------ ------ ------ ------- -------- -------- ---- -------- ------- ------ -------- ------ -------- ----- ------ -------- ------ -------- -------- -------- -------- ------
    20140426094659 2014-04-26T09:46:59.000Z 2014-04-26T09:53:41.300Z 1564.0199 1568.0461        3    None            None    None   None  None          None    None    None    None    None    None    None    None    None    None    None    None    None    None    None    None    None   None              None 2014-04-26T00:00:00.000Z None    None     None     None    None                                None     None  None                     None    None  None  None  None  None  None  None  None  None None             5    105    471   1010   None    None     None     None None     None    None   None            None     None  None   None     None   None   GREGOR     None     None     1566      0
    20140426095643 2014-04-26T09:56:43.000Z 2014-04-26T10:09:50.300Z 1564.0191 1568.0461        3    None            None    None   None  None          None    None    None    None    None    None    None    None    None    None    None    None    None    None    None    None    None   None              None 2014-04-26T00:00:00.000Z None    None     None     None    None                                None     None  None                     None    None  None  None  None  None  None  None  None  None None             5    200    471   1010   None    None     None     None None     None    None   None            None     None  None   None     None   None   GREGOR     None     None     1566      0
    <BLANKLINE>
    <BLANKLINE>
    """

    @property
    def info_url(self):
        return 'https://solarnet2.oma.be'



    def search(self, *query):
        """
        Searches the SVO for datasets matching the query.

        Parameters
        ----------
        *query : `~sunpy.net.attrs`
            Attributes to search for datasets.

        Returns
        -------
        `~sunpy.net.base_client.QueryResponseTable`
            A table containing the search results.

        Examples
        --------
        >>> from sunpy.net import Fido, attrs as a
        >>> query = [a.solarnet.Dataset.swap_level_1, a.solarnet.Limit(1)]
        >>> search_results = Fido.search(*query)   # doctest: +REMOTE_DATA
        >>> search_results  # doctest: +REMOTE_DATA
        <sunpy.net.fido_factory.UnifiedResponse object at ...>
        Results from 1 Provider:
        <BLANKLINE>
        1 Results from the SOLARNETClient:
        Source: https://solarnet2.oma.be
        <BLANKLINE>
            Oid               Date_beg                 Date_end         Wavemin Wavemax Artefx Bscale Bunit Bzero  Cap_mode   Cd1_1   Cd1_2 Cd2_1   Cd2_2     Cdelt1    Cdelt2  Compress      Creator       Crota1 Crota2 Crpix1 Crpix2 Crval1 Crval2  Ctype1   Ctype2  Cunit1 Cunit2 Datamax Datamin           Date                   Date_obs         Detector    Dsun_obs    Dtplar1 Dtplar2    Eacqtime   Exptime            Filename           Filter Firstcol Firstrow    Geod_alt      Geod_lat      Geod_lon      Gsex_obs       Gsey_obs      Gsez_obs   Hasblack Hasoffst Hasstdby    Heex_obs       Heey_obs      Heez_obs       Hgln_obs        Hglt_obs   Instrume Is_proc Lang_rot Last_col Last_row Led_pow Led_sel Level Lonpole    Los_alt    Lzwdecor Naxis Naxis1 Naxis2 Nprescr Npreslzw  Object   Obs_mode   Origin P2_roll P2_x0 P2_y0 Pav_rot0 Pav_rot1 Pga_gain Pga_offs  On Readrdiv Rebin Recbias Recnum Recoding    Rsun_arc      Sacqtime   Sizcompi    Solar_ep   Swavint Swxcen Swycen Telescop    Temp1det       Temp2det       Tempdark       Trantime   Trapelec Trapprot          Ttemp1                  Ttemp2         Wavelnth          Wcsname           F_size
        -------------- ------------------------ ------------------------ ------- ------- ------ ------ ----- ------ -------- --------- ----- ----- --------- --------- --------- -------- ------------------ ------ ------ ------ ------ ------ ------ -------- -------- ------ ------ ------- ------- ------------------------ ------------------------ -------- -------------- ------- ------- ------------- ------- ----------------------------- ------ -------- -------- ------------- ------------- ------------- -------------- ------------- ------------- -------- -------- -------- -------------- ------------- ------------- ---------------- ------------- -------- ------- -------- -------- -------- ------- ------- ----- ------- ------------- -------- ----- ------ ------ ------- -------- ------- ------------ ------ ------- ----- ----- -------- -------- -------- -------- --- -------- ----- ------- ------ -------- ------------- ------------- -------- ------------- ------- ------ ------ -------- -------------- -------------- -------------- ------------- -------- -------- ----------------------- ----------------------- -------- ------------------------- -------
        20091120082529 2009-11-20T08:25:29.027Z 2009-11-20T08:25:30.027Z    17.4    17.4    off 0.0625  DN/s 2048.0       DS 3.1646941   0.0   0.0 3.1646941 3.1646941 3.1646941      off P2SW_PREP.PRO v1.5    0.0    0.0  512.5  512.5    0.0    0.0 HPLN-TAN HPLT-TAN arcsec arcsec 50.3125     0.0 2017-11-16T15:30:05.000Z 2009-11-20T08:25:29.027Z     SWAP 147819620849.0  2000.0   177.0 25691981838.0     1.0 swap_lv1_20091120_082529.fits     Al        1        1 723227.647021 49.8713348902 153.634430938 -2841004.71812 3290864.59669 5599236.25546        4       11        0 147819620706.0 -3290864.5967 5599236.25545 -0.0010059134686 2.20683198122     SWAP       0      0.0     1024     1024     off       a     1   180.0 126023.112807      off     2   1024   1024       0        0 Sun EUV Sun centered    ROB    None  None  None      0.0      0.0        1       59   0        0   off       0      0      off 971.178544869 25691900864.0        0 353.089769721 21.4853 353.01 447.68   PROBA2 -2.88999023437 -4.00998535156 -4.00998535156 25691981838.0      0.0      0.0 2009-11-20T08:09:05.000 2009-11-20T08:28:58.000      174 Helioprojective-cartesian 2113920
        <BLANKLINE>
        <BLANKLINE>

        """
        results = []
        query = and_(*query)
        block = walker.create(query)[0]

        if "datasets" in block:
            url = _BASE_URL.format(block["datasets"])
            block.pop("datasets")

        try:
            req = urllib.request.urlopen(url + '?' + urllib.parse.urlencode(block))
            data = json.loads(req.read().decode())["objects"]
        except Exception as e:
            raise Exception(f"Failed to fetch data from {url}: {str(e)}")

        hide_keys = ["Fits_header", "Data_location", "Resource_uri", "url", "File_raw", "File_tar", "File_tmr"]

            # just to avoid empty tags
        hide_keys.append("Tags") if "tag" not in block else hide_keys

        for i in data:
                    # Filter out keys that should be hidden or not displayed
                filtered_data = {key.capitalize(): value for key, value in i.items()}

                    # Add specific fields from data_location
                filtered_data.update({
                        "F_size": i["data_location"]["file_size"],
                        "url": i["data_location"]["file_url"],
                    })

                results.append(filtered_data)



        qrt = QueryResponseTable(results, client=self)
        qrt.hide_keys = hide_keys
        return qrt

    @staticmethod
    def load_solarnet_values():
        """
        Load Solarnet dataset values from a local JSON file.

        Returns
        -------
        attrs : `dict`
            Dictionary of Solarnet dataset values.
        """
        data_sets = Path(__file__).parent / "data" / "datasets.json"
        tags = Path(__file__).parent / "data" / "tags.json"

        with data_sets.open() as data_values:
            data = json.load(data_values)

        with tags.open() as tags_values:
            tags = json.load(tags_values)

        data_set = {Dataset: list(data.items())}
        target = {Target: [("AR", "Active Region"),("CH","Coronal Hole"),("FS","Flare"),("QR","Quiet Region")]}
        tags = {Tags:list(tags.items())}

        attrs = data_set | target | tags
        return attrs

    def create_parse_solarnet_values():
        """
        Gets the names of available datasets, tags, and targets from SolarNet and saves them into separate JSON files.
        """
        dir_path = os.path.dirname(os.path.realpath(__file__))
        data_dir = os.path.join(dir_path, 'data')
        os.makedirs(data_dir, exist_ok=True)

        # Datasets
        url_dataset = _BASE_URL.format("dataset")
        response = requests.get(url_dataset, params={"limit": 100})
        data = response.json()
        values = {}
        for obj in data.get("objects", []):
            name = obj["name"].replace(" ", "_").lower()
            description = obj.get("description") or obj.get("telescope", {}).get("description", "")
            values[name] = description.split(". ")[0].split(", ")[0] if description else ""
        with open(os.path.join(data_dir, 'datasets.json'), 'w') as f:
            json.dump(dict(sorted(values.items())), f, indent=2)

        # Tags
        url_tags = _BASE_URL.format("tag")
        response = requests.get(url_tags)
        data_tags = response.json()
        tag_values = {}
        for obj in data_tags.get("objects", []):
            tag_name = obj["name"].replace(" ", "_").lower()
            tag_values[tag_name] = tag_name
        with open(os.path.join(data_dir, 'tags.json'), 'w') as f:
            json.dump(dict(sorted(tag_values.items())), f, indent=2)

    def fetch(self, query_results, path, downloader, **kwargs):
        for row in query_results:
            fname = row['url'].split('/')[-1]
            filepath = str(path).format(file=fname)
            max_splits = kwargs.get('max_splits', 5)
            downloader.enqueue_file(row['url'], filename=filepath, max_splits=max_splits)

    @classmethod
    def register_values(cls):
        return cls.load_solarnet_values()

    @classmethod
    def _attrs_module(cls):
        return 'solarnet', 'sunpy.net.solarnet.attrs'

    @classmethod
    def _can_handle_query(cls, *query):
        return any(isinstance(q, a.solarnet.Dataset) for q in query)
