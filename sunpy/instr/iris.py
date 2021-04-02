from sunpy.extern.sunkit_instruments.iris import SJI_to_sequence

__all__ = ['SJI_to_sequence']

# Trick the docs into thinking these functions are defined in here.
for _a in (SJI_to_sequence,):
    _a.__module__ = __name__
