from __future__ import absolute_import

__all__ = ['OrderedDict']

from astropy.utils.compat.odict import OrderedDict

def __meta_repr__(self):
    '''This function overrides the default __repr__ function for
    OrderedDict to make the output more human-readable. This is
    particularly relevant for displaying lightcurve.meta information.'''
    totalstring=''
    for l,d in self.items():
        out = (str(l) + ':\t' + str(d) + '\n').expandtabs(20)
        totalstring=totalstring + out
            
    return (totalstring)


OrderedDict.__repr__ = __meta_repr__
