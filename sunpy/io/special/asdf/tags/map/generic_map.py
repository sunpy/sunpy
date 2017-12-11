# -*- coding: utf-8 -*-
import numpy as np

import astropy.units as u

import sunpy.map
from ...types import SunPyType

__all__ = ['GenericMapType']


class GenericMapType(SunPyType):
    name = "map/generic_map"
    types = ['sunpy.map.GenericMap']
    requires = ['sunpy']
    version = "1.0.0"

    @classmethod
    def from_tree(cls, node, ctx):
        # Use the factory here to get the correct subclass back
        out_map = sunpy.map.Map(node['data'], node['header'])
        out_map.shift(*node['shift'])
        return out_map

    @classmethod
    def to_tree(cls, smap, ctx):
        node = {}
        node['data'] = np.asarray(smap.data)
        node['meta'] = dict(smap.meta)
        node['shift'] = u.Quantity(smap.shifted_value)

        # TODO: Save some or all of plot_settings
        # node['plot_settings'] = smap.plot_settings

        return node
