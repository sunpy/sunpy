import json
import pathlib
from concurrent.futures import ThreadPoolExecutor

import requests

from astropy.table import Table

from .cdaweb import _CDAS_BASEURL, _CDAS_HEADERS, _DATAVIEW

__all__ = ['get_observatory_groups', 'get_datasets']


def get_observatory_groups():
    """
    Get a list of observatory IDs for each observatory in CDAWeb.

    An observatory group is typically a single mission, which can contain
    multiple observatories, e.g. for the STEREO observatory group there are two
    observatories, STEREO-A and STEREO-B.

    Returns
    -------
    `astropy.table.Table`

    Examples
    --------
    >>> from sunpy.net.cdaweb import get_observatory_groups
    >>>
    >>> groups = get_observatory_groups() # doctest: +REMOTE_DATA
    >>> groups['Group'] # doctest: +REMOTE_DATA
        <Column name='Group' dtype='str55' length=...>
                        ACE
                        AIM
                      AMPTE
        ...
                    Voyager
                       Wind
    >>> groups.loc['STEREO'] # doctest: +REMOTE_DATA
    <Row index=...>
    Group                                  Observatories
    str55                                      str...
    ------ -----------------------------------------------------------------------------
    STEREO 'Ahead', 'Behind', 'STA', 'STB', 'STEREO', 'STEREOA', 'STEREOB', 'sta', 'stb'
    """
    # Get a list of files for a given dataset between start and end times
    url = '/'.join([
        _CDAS_BASEURL,
        'dataviews', _DATAVIEW,
        'observatoryGroups'
    ])
    response = requests.get(url, headers=_CDAS_HEADERS)
    obs_groups = response.json()

    names = [obs['Name'] for obs in obs_groups['ObservatoryGroupDescription']]
    obs_ids = [obs['ObservatoryId'] for obs in obs_groups['ObservatoryGroupDescription']]
    # Join all IDs into a single string
    obs_ids = ["'" + "', '".join(id) + "'" for id in obs_ids]

    t = Table([names, obs_ids], names=['Group', 'Observatories'])
    t.add_index('Group')
    return t


def get_datasets(observatory):
    """
    Get a list of datasets for a given observatory.

    Parameters
    ----------
    observatory : `str`
        Observatory name.

    Returns
    -------
    `astropy.table.Table`

    Examples
    --------
    >>> from sunpy.net.cdaweb import get_datasets
    >>>
    >>> datasets = get_datasets('STEREOB') # doctest: +REMOTE_DATA
    >>> datasets['Id'] # doctest: +REMOTE_DATA
    <Column name='Id' dtype='str17' length=5>
        STB_LB_IMPACT
    STB_L1_IMPACT_HKP
           STB_L1_HET
      STB_L2_SWEA_PAD
     STB_L1_SWEA_SPEC
    >>> datasets.loc['STB_L1_SWEA_SPEC']['Label'] # doctest: +REMOTE_DATA
    np.str_('STEREO Behind IMPACT/SWEA Spectra - J. Luhmann (UCB/SSL)')
    >>> datasets.loc['STB_L1_SWEA_SPEC'][['Start', 'End']] # doctest: +REMOTE_DATA
    <Row index=4>
             Start                     End
             str24                    str24
    ------------------------ ------------------------
    2012-12-01T00:00:03.000Z 2013-12-31T23:59:41.000Z
    """
    # Get a list of files for a given dataset between start and end times
    url = '/'.join([
        _CDAS_BASEURL,
        'dataviews', _DATAVIEW,
        'datasets'
    ])
    url = f'{url}?observatory={observatory}'
    response = requests.get(url, headers=_CDAS_HEADERS)
    datasets = response.json()['DatasetDescription']

    ids = [dataset['Id'] for dataset in datasets]
    instruments = [', '.join(dataset['Instrument']) for dataset in datasets]
    labels = [dataset['Label'] for dataset in datasets]
    stimes = [dataset['TimeInterval']['Start'] for dataset in datasets]
    etimes = [dataset['TimeInterval']['End'] for dataset in datasets]

    t = Table([ids, instruments, labels, stimes, etimes],
              names=['Id', 'Instruments', 'Label', 'Start', 'End'])
    t.add_index('Id')
    return t


def _update_cdaweb_dataset_data():
    all_obs = get_observatory_groups()
    url = '/'.join([
        _CDAS_BASEURL,
        'dataviews', _DATAVIEW,
        'datasets'
    ])
    # Mapping from dataset ID to description
    all_datasets = {}
    # Number of parallel threads we spawn
    N = 3

    def _fetch_cdaweb_dataset(group, url=url):
        print(f'🛰 Getting datasets for {group}')
        u = url + f'?observatoryGroup={group}'
        res = requests.get(u, headers=_CDAS_HEADERS)
        datasets = res.json()['DatasetDescription']
        dataset_ids = {ds['Id']: ds['Label'] for ds in datasets}
        all_datasets.update(dataset_ids)

    with ThreadPoolExecutor(max_workers=N) as executor:
        # Submit each URL to the thread pool
        futures = [executor.submit(_fetch_cdaweb_dataset, group) for group in all_obs['Group']]
        # Wait for all tasks to complete
        for future in futures:
            future.result()

    attr_file = pathlib.Path(__file__).parent / 'data' / 'attrs.json'
    with open(attr_file, 'w') as attrs_file:
        json.dump(dict(sorted(all_datasets.items())), attrs_file, indent=2)
