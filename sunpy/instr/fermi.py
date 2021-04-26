from sunpy.extern.sunkit_instruments.fermi import (
    download_weekly_pointing_file,
    get_detector_sun_angles_for_date,
    get_detector_sun_angles_for_time,
    met_to_utc,
    plot_detector_sun_angles,
)

__all__ = ['download_weekly_pointing_file', 'get_detector_sun_angles_for_time',
           'get_detector_sun_angles_for_date', 'plot_detector_sun_angles',
           'met_to_utc']

# Trick the docs into thinking these functions are defined in here.
for _a in (download_weekly_pointing_file,
           get_detector_sun_angles_for_date,
           get_detector_sun_angles_for_time,
           met_to_utc,
           plot_detector_sun_angles,):
    _a.__module__ = __name__
