import numpy as np

import astropy.units as u
from asdf.yamlutil import custom_tree_to_tagged_tree

import sunpy.map
from sunpy.io.special.asdf.types import SunPyType

__all__ = ['GenericMapType']


class GenericMapType(SunPyType):
    name = "map/generic_map"
    types = ['sunpy.map.GenericMap']
    requires = ['sunpy']
    version = "1.0.0"

    @classmethod
    def from_tree(cls, node, ctx):
        # Use the factory here to get the correct subclass back
        out_map = sunpy.map.Map(np.asarray(node['data']), node['meta'])
        out_map.shift(*node['shift'])

        out_map.mask = node.get('mask')
        out_map.uncertainty = node.get('uncertainty')
        out_map._unit = node.get('unit')

        return out_map

    @classmethod
    def to_tree(cls, smap, ctx):
        node = {}
        node['data'] = np.asarray(smap.data)
        node['meta'] = dict(smap.meta)
        node['shift'] = u.Quantity(smap.shifted_value)
        node['mask'] = smap.mask
        node['uncertainty'] = smap.uncertainty
        node['unit'] = smap.unit

        # TODO: Save some or all of plot_settings
        # node['plot_settings'] = smap.plot_settings

        return custom_tree_to_tagged_tree(node, ctx)

    @classmethod
    def assert_equal(cls, old, new):
        """
        This method is used by asdf to test that to_tree > from_tree gives an
        equivalent object.
        """
        np.testing.assert_allclose(old.data, new.data)

        # Test the meta by force!
        for ok, ov in old.meta.items():
            assert ok in new.meta
            assert new.meta[ok] == ov

        assert u.allclose(old.shifted_value, new.shifted_value)
        if old.mask is not None and new.mask is not None:
            np.testing.assert_allclose(old.mask, new.mask)
        assert old.unit == new.unit
