from sunpy.net import attrs as a
from sunpy.net.dataretriever import GenericClient
from sunpy.net.dataretriever.attrs.adapt import (
    ADAPTDataAssimilation,
    ADAPTEvolutionMode,
    ADAPTFileType,
    ADAPTHelioData,
    ADAPTInputSource,
    ADAPTLonType,
    ADAPTMagData,
    ADAPTRealizations,
    ADAPTResolution,
    ADAPTVersionMonth,
    ADAPTVersionYear,
)

__all__ = ['ADAPTClient']


class ADAPTClient(GenericClient):
    """
    Provides access to the ADvanced Adaptive Prediction Technique (ADAPT) products
    of the National Solar Observatory (NSO).

    `Searches data hosted by the NSO <https://gong.nso.edu/adapt/maps/gong/>`__

    Examples
    --------
    >>> import astropy.units as u

    >>> from sunpy.net import Fido, attrs as a
    >>> from sunpy.coordinates.sun import carrington_rotation_time

    >>> # Define the Carrington Rotation Number and the number of frames
    >>> CR = 2193
    >>> frames = 10
    >>> date_start = carrington_rotation_time(CR)
    >>> date_end = date_start + frames*(3*1.9999999 * u.hour)
    >>> longitude_type = '0'

    >>> Fido.search(a.Time(date_start, date_end), a.Instrument('adapt'), a.adapt.ADAPTLonType(longitude_type))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    10 Results from the ADAPTClient:
    <BLANKLINE>
           Start Time               End Time        Instrument Provider Source ... ADAPTMagData days_since_last_obs hours_since_last_obs minutes_since_last_obs seconds_since_last_obs
    ----------------------- ----------------------- ---------- -------- ------ ... ------------ ------------------- -------------------- ---------------------- ----------------------
    2017-07-20 08:00:00.000 2017-07-20 08:00:59.999      ADAPT      NSO   GONG ...            1                   0                    1                     56                      0
    2017-07-20 14:00:00.000 2017-07-20 14:00:59.999      ADAPT      NSO   GONG ...            1                   0                    1                     56                      0
    2017-07-20 20:00:00.000 2017-07-20 20:00:59.999      ADAPT      NSO   GONG ...            1                   0                    1                     56                      0
    2017-07-21 02:00:00.000 2017-07-21 02:00:59.999      ADAPT      NSO   GONG ...            1                   0                    1                     56                      0
    2017-07-21 08:00:00.000 2017-07-21 08:00:59.999      ADAPT      NSO   GONG ...            1                   0                    1                     56                      0
    2017-07-21 14:00:00.000 2017-07-21 14:00:59.999      ADAPT      NSO   GONG ...            1                   0                    1                     56                      0
    2017-07-21 20:00:00.000 2017-07-21 20:00:59.999      ADAPT      NSO   GONG ...            1                   0                    1                     56                      0
    2017-07-22 02:00:00.000 2017-07-22 02:00:59.999      ADAPT      NSO   GONG ...            1                   0                    1                     56                      0
    2017-07-22 08:00:00.000 2017-07-22 08:00:59.999      ADAPT      NSO   GONG ...            1                   0                    4                     36                      0
    2017-07-22 14:00:00.000 2017-07-22 14:00:59.999      ADAPT      NSO   GONG ...            1                   0                    1                     56                      0
    <BLANKLINE>
    <BLANKLINE>

    References
    ----------
    `Names and possible attrs values are available <https://gong.nso.edu/adapt/maps/adapt_filename_notes.txt>`__.
    """
    pattern = r'https://gong.nso.edu/adapt/maps/gong/{{year:4d}}/adapt{{ADAPTFileType:1d}}{{ADAPTLonType:1d}}{{ADAPTInputSource:1d}}{{ADAPTDataAssimilation:1d}}{{ADAPTResolution:1d}}' + \
    '_{{ADAPTVersionYear:2d}}{{ADAPTVersionMonth:1l}}{{ADAPTRealizations:3d}}_{{year:4d}}{{month:2d}}{{day:2d}}{{hour:2d}}{{minute:2d}}' + \
    '_{{ADAPTEvolutionMode:1l}}{{days_since_last_obs:2d}}{{hours_since_last_obs:2d}}{{minutes_since_last_obs:2d}}{{seconds_since_last_obs:2d}}{{ADAPTHelioData:1l}}{{ADAPTMagData:1d}}.fts.gz'

    @classmethod
    def _attrs_module(cls):
        return 'adapt', 'sunpy.net.dataretriever.attrs.adapt'

    @classmethod
    def register_values(cls):
        adict = {
            a.Instrument: [('ADAPT', 'ADvanced Adaptive Prediction Technique.')],
            a.Provider: [('NSO', 'National Solar Observatory.')],
            a.Source: [('GONG', 'Global Oscillation Network Group.')],
            ADAPTFileType: [('4', 'Public')],
            ADAPTLonType: [('0', 'Carrington Fixed'), ('1', 'Central Meridian'), ('2', 'East Limb')],
            ADAPTInputSource: [('0', 'All'), ('1', 'KPVT'), ('2', 'VSM'), ('3', 'GONG'), ('4', 'HMI'), ('5', 'MDI'), ('6', 'MWO')],
            ADAPTDataAssimilation: [('0', 'WH'), ('1', 'enLS'), ('2', 'enkf'), ('3', 'enLAKF')],
            ADAPTResolution: [('1', '1.0 deg'), ('2', '0.2 deg')],
            ADAPTVersionYear: [(str(i), f"Code version year -> {2000 + i}") for i in range(1, 20)],
            ADAPTRealizations: [(str(i), f"Number of Realizations -> {i}") for i in range(1, 20)],
            ADAPTVersionMonth: [(chr(i+96), f"Code version month -> {i}") for i in range(1, 13)],
            ADAPTEvolutionMode: [('a', 'Data assimilation step'), ('i', 'Intermediate step'), ('f', 'Forecast step')],
            ADAPTHelioData: [('n', 'Not added or no data'), ('f', 'Far-side'), ('e', 'Emergence'), ('b', 'Both emergence & far-side')],
            ADAPTMagData: [('0', 'Not added or no data'), ('1', 'Mag-los'), ('2', 'Mag-vector'), ('3', 'Mag- both los & vector'),
                            ('4', 'Mag- polar avg obs'), ('5', 'Mag- los & polar'), ('6', 'Mag- vector & polar'), ('7', 'Mag- both los and vector & polar')]
        }
        return adict
