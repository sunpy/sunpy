from __future__ import absolute_import

__all__ = ['OrderedDict']

from astropy.utils.compat.odict import OrderedDict

def __meta_str__(self):
    '''This function overrides the default __str__ function for
    OrderedDict to make the output more human-readable. This is
    particularly relevant for displaying lightcurve.meta information.'''
    out=[]
    for l,d in self.iteritems():
        out.append( ("""%s: \t %s""" % (l,d)).expandtabs(20))
    totalstring = '\n'.join(out)
            
    return totalstring


OrderedDict.__str__ = __meta_str__
