from sunpy.net import attrs
from sunpy.net.fido_factory import Fido
from sunpy.net.dataretriever import GenericClient
from sunpy.net.dataretriever.attrs.adapt import *
import sunpy.net.dataretriever.attrs as atrs
from sunpy.net.attr import SimpleAttr


__all__ = ['ADAPTClient']


class ADAPTClient(GenericClient):
    """
    Provides access to the ADvanced Adaptive Prediction Technique (ADAPT) products of the National Solar Observatory.

        Searches data hosted by the `National Solar Observatory <https://gong.nso.edu/adapt/maps/gong/>`__

        Examples
        --------
        >>> from sunpy.net import Fido, attrs as a
        >>> get_date, get_date_end = adapt_source.carrington_time(CR, frames)
        >>> LngType = '0' # 0 is carrington, 1 is central meridian
        >>> print (f"Searching for images from {get_date} to {get_date_end}\n")
        >>> res = Fido.search(attrs.Instrument('adapt'), attrs.Time(get_date, get_date_end), ADAPTLngType(LngType))
        >>> print(res)
        >>> ret =  Fido.fetch(res, path=make_dir())


    Searching for images from 2017-07-20T06:03:05 to 2017-07-20T12:03:05

    Results from 1 Provider:

    1 Results from the ADAPTClient:

        Start Time               End Time        ... minutes_since_last_obs seconds_since_last_obs
    ----------------------- ----------------------- ... ---------------------- ----------------------
    2017-07-20 08:00:00.000 2017-07-20 08:00:59.999 ...                     56                      0


    Fetching search_result:
    Files Downloaded: 100%|███████████████████████████████████████████| 1/1 [00:00<00:00,  7.65file/s]
    Success!
    """
    baseurl = r'https://gong.nso.edu/adapt/maps/gong/%Y/adapt(\d){5}_(\d){2}(\w){1}(\d){3}_(\d){12}_(\w){1}(\d){8}(\w){1}(\d){1}\.fts\.gz'
    pattern = '{}adapt{ADAPTFileType:1d}{ADAPTLngType:1d}{ADAPTInputSource:1d}{ADAPTDataAssimilation:1d}{ADAPTResolution:1d}_{ADAPTVersionYear:2d}{ADAPTVersionMonth:1l}{ADAPTRealizations:3d}_{year:4d}{month:2d}{day:2d}{hour:2d}{minute:2d}_{ADAPTEvolutionMode:1l}{days_since_last_obs:2d}{hours_since_last_obs:2d}{minutes_since_last_obs:2d}{seconds_since_last_obs:2d}{ADAPTHelioData:1l}{ADAPTMagData:1d}.fts.gz'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs as a

        adict ={attrs.Instrument: [('ADAPT', 'ADvanced Adaptive Prediction Technique.')],
                attrs.Source: [('NSO', 'National Solar Observatory.')],
                attrs.Provider: [('GONG', 'Global Oscillation Network Group.')],
                ADAPTFileType: [('4', 'Public')],
                ADAPTLngType: [('0', 'Carrington Fixed'), ('1', 'Central Meridian'), ('2', 'East Limb')],
                ADAPTInputSource: [('0', 'All'), ('1', 'KPVT'), ('2', 'VSM'), ('3', 'GONG'), ('4', 'HMI'), ('5', 'MDI'), ('6', 'MWO')],
                ADAPTDataAssimilation: [('0', 'WH'), ('1', 'enLS'), ('2', 'enkf'), ('3', 'enLAKF')],
                ADAPTResolution: [('1', '1.0 deg'), ('2', '0.2 deg')],
                ADAPTVersionYear: [(str(i), f"Code version year -> {2000 + i}") for i in range(1, 20)],
                ADAPTRealizations: [(str(i), f"Number of Realizations -> {i}") for i in range(1, 20)],
                ADAPTVersionMonth: [(chr(i+96), f"Code version month -> {i}") for i in range(1, 13)],
                ADAPTEvolutionMode: [('a', 'Data assimilation step'), ('i', 'Intermediate step'), ('f', 'Forecast step')],
                ADAPTHelioData: [('n', 'Not added or no data'), ('f', 'Far-side'), ('e', 'Emergence'), ('b', 'Both emergence & far-side')],
                ADAPTMagData: [('0', 'Not added or no data'), ('1', 'Mag-los'), ('2', 'Mag-vector'), ('3', 'Mag- both los & vector'), ('4', 'Mag- polar avg obs'), ('5', 'Mag- los & polar'), ('6', 'Mag- vector & polar'), ('7', 'Mag- both los and vector & polar')]
                }
        return adict

    @classmethod
    def _can_handle_query(cls, *query):
        required = {attrs.Instrument, attrs.Time}

        optional = {ADAPTFileType, ADAPTLngType, ADAPTInputSource, ADAPTDataAssimilation,
                    ADAPTResolution, ADAPTVersionYear, ADAPTVersionMonth, ADAPTEvolutionMode,
                    ADAPTHelioData, ADAPTMagData}

        all_attrs = {type(x) for x in query}
        # print(all_attrs)
        return required.issubset(all_attrs) and all_attrs.issubset(required.union(optional))

def carrington_time(CR=2193, frames=1):
    """Get the start date for the start of a carrington rotation,
    and a duration that will retrieve a given number of frames

    Returns:
        _type_: _description_
    """
    import astropy.units as u
    import sunpy.coordinates
    date = sunpy.coordinates.sun.carrington_rotation_time(CR)
    date_end = date + frames*(3*1.9999999 * u.hour)

    # Format the Search Dates
    tstring = r"%Y-%m-%dT%H:%M:%S"
    get_date    =date.strftime(tstring)
    get_date_end=date_end.strftime(tstring)
    return get_date, get_date_end


if __name__ == "__main__":
    from sunpy.net.dataretriever.sources.tests.test_adapt import test_client
    test_client()
