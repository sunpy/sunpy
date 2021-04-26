from sunpy.extern.sunkit_instruments.lyra import (
    get_lytaf_event_types,
    get_lytaf_events,
    remove_lytaf_events_from_timeseries,
    split_series_using_lytaf,
)

__all__ = ['remove_lytaf_events_from_timeseries',
           'get_lytaf_events',
           'get_lytaf_event_types',
           'split_series_using_lytaf']

# Trick the docs into thinking these functions are defined in here.
for _a in (get_lytaf_event_types,
           get_lytaf_events,
           remove_lytaf_events_from_timeseries,
           split_series_using_lytaf):
    _a.__module__ = __name__
