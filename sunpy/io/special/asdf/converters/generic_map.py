import numpy as np

from asdf.extension import Converter

__all__ = ["GenericMapConverter"]


class GenericMapConverter(Converter):
    tags = ["tag:sunpy.org:sunpy/map/generic_map-*"]

    @property
    def types(self):
        from sunpy.map.mapbase import GenericMap
        return [GenericMap] + list(GenericMap._registry.keys())

    def select_tag(self, obj, tags, ctx):
        # Sort the tags in reverse alphabetical order and pick the first (i.e.
        # the one with the highest version). This assumes that all the tags for
        # this converter are named the same other than the version number.
        tags = list(sorted(tags, reverse=True))
        return tags[0]

    def from_yaml_tree(self, node, tag, ctx):
        import astropy.units as u

        import sunpy.map

        # Use the factory here to get the correct subclass back
        out_map = sunpy.map.Map(np.asanyarray(node["data"]), node["meta"])
        out_map.mask = node.get("mask")
        out_map.uncertainty = node.get("uncertainty")
        out_map._unit = node.get("unit")

        # Old files can have a shift that needs applying
        if 'generic_map-1.0.0' in tag and not np.allclose(node['shift'], 0 * u.deg):
            out_map = out_map.shift_reference_coord(*node['shift'])

        return out_map

    def to_yaml_tree(self, smap, tag, ctx):
        node = {}
        node["data"] = np.asarray(smap.data)
        node["meta"] = dict(smap.meta)
        if smap.mask is not None:
            node["mask"] = smap.mask
        if smap.uncertainty is not None:
            node["uncertainty"] = smap.uncertainty
        if smap.unit is not None:
            node["unit"] = smap.unit

        # TODO: Save some or all of plot_settings
        # node["plot_settings"] = smap.plot_settings

        return node
