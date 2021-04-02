from sunpy.extern.sunkit_instruments.goes_xrs import (
    calculate_radiative_loss_rate,
    calculate_temperature_em,
    calculate_xray_luminosity,
    flareclass_to_flux,
    flux_to_flareclass,
    get_goes_event_list,
)

__all__ = ['get_goes_event_list', 'calculate_temperature_em',
           'calculate_radiative_loss_rate', 'calculate_xray_luminosity', 'flux_to_flareclass',
           'flareclass_to_flux']

# Trick the docs into thinking these functions are defined in here.
for _a in (calculate_radiative_loss_rate,
           calculate_temperature_em,
           calculate_xray_luminosity,
           flareclass_to_flux,
           flux_to_flareclass,
           get_goes_event_list):
    _a.__module__ = __name__
