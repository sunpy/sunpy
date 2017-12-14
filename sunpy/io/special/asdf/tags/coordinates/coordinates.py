# -*- coding: utf-8 -*-

from asdf.yamlutil import custom_tree_to_tagged_tree

import astropy.units as u
from astropy.coordinates import Longitude, Latitude, Angle
from astropy.io.misc.asdf.tags.unit.quantity import QuantityType

from sunpy.coordinates import HeliographicStonyhurst, HeliographicCarrington
from sunpy.io.special.asdf.types import SunPyType


__all__ = ['HGCoordType']


class HGCoordType(SunPyType):
    name = "coordinates/hgs_coord"
    types = ['sunpy.coordinates.HeliographicStonyhurst']
    requires = ['sunpy']
    version = "1.0.0"

    @classmethod
    def to_tree(cls, frame, ctx):
        node = {}

        node['name'] = frame.name
        wrap_angle = u.Quantity(frame.Tx.wrap_angle)
        node['lon'] = {
            'value': frame.lon.value,
            'unit': frame.lon.unit.to_string(),
            'wrap_angle': custom_tree_to_tagged_tree(wrap_angle, ctx)
        }
        node['lat'] = custom_tree_to_tagged_tree(u.Quantity(frame.lat))
        node['distance'] = custom_tree_to_tagged_tree(u.Quantity(frame.distance))
        node['obstime'] = custom_tree_to_tagged_tree(frame.obstime, ctx)

        return node

    @classmethod
    def from_tree(cls, node, ctx):
        angle = Angle(QuantityType.from_tree(node['lat']['wrap_angle'], ctx))
        wrap_angle = Angle(angle)
        name = node['name']
        lon = Longitude(node['lon']['value'],
                        unit=node['lon']['unit'],
                        wrap_angle=wrap_angle)
        lat = Latitude(node['lat'])
        distance = u.Quantity(node['distance'])
        obstime = node['obstime']

        if name == 'heliographic_stonyhurst':
            return HeliographicStonyhurst(lon=lon, lat=lat, distance=distance, obstime=obstime)
        elif name == 'heliographic_carrington':
            return HeliographicCarrington(lon=lon, lat=lat, distance=distance, obstime=obstime)
        else:
            raise ValueError("Unknown Heliographic coordinate name: {}".format(name))
