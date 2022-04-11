import numpy as np

import astropy.units as u
from asdf.extension import Converter

__all__ = ["GenericMapConverter"]


class GenericMapConverter(Converter):
    tags = ["tag:sunpy.org:sunpy/map/generic_map-*"]

    @property
    def types(self):
        from sunpy.map.mapbase import GenericMap
        return list(GenericMap._registry.keys())

    def from_yaml_tree(self, node, tag, ctx):
        import sunpy.map

        # Use the factory here to get the correct subclass back
        out_map = sunpy.map.Map(np.asarray(node["data"]), node["meta"])
        out_map = out_map.shift(*node["shift"])

        out_map.mask = node.get("mask")
        out_map.uncertainty = node.get("uncertainty")
        out_map._unit = node.get("unit")
        return out_map

    def to_yaml_tree(self, smap, tag, ctx):
        node = {}
        node["data"] = np.asarray(smap.data)
        node["meta"] = dict(smap.meta)
        node["shift"] = u.Quantity(smap.shifted_value)
        if smap.mask is not None:
            node["mask"] = smap.mask
        if smap.uncertainty is not None:
            node["uncertainty"] = smap.uncertainty
        if smap.unit is not None:
            node["unit"] = smap.unit

        # TODO: Save some or all of plot_settings
        # node["plot_settings"] = smap.plot_settings

        return node
