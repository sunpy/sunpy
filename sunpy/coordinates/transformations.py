# -*- coding: utf-8 -*-
# NumPy import
import numpy as np

# Astropy imports
from astropy import units as u
from astropy.coordinates.representation import CartesianRepresentation
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import FunctionTransform

# SunPy imports
from sunpy import sun

from .representation import SphericalWrap180Representation
from .frames import (HelioGraphicStonyhurst, HelioGraphicCarrington,
                     HelioCentric, HelioProjective)


def _carrington_offset(dateobs):
    """
    Calculate the HG Longitude offest based on a time
    """
    if dateobs is None:
        raise ValueError("To perform this transformation the coordinate Frame needs a dateobs Attribute")
    return sun.heliographic_solar_center(dateobs)[0]


#==============================================================================
#------------------------- Transformation Framework ---------------------------
#==============================================================================


@frame_transform_graph.transform(FunctionTransform, HelioGraphicStonyhurst, HelioGraphicCarrington)
def hgs_to_hgc(hgscoord, hgcframe):
    c_lon = hgscoord.spherical.lon + _carrington_offset(hgscoord.dateobs).to(u.deg)
    representation = SphericalWrap180Representation(c_lon, hgscoord.lat, hgscoord.rad)
    return hgcframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, HelioGraphicCarrington, HelioGraphicStonyhurst)
def hgc_to_hgs(hgccoord, hgsframe):
    s_lon = hgccoord.spherical.lon - _carrington_offset(hgccoord.dateobs).to(u.deg)
    representation = SphericalWrap180Representation(s_lon, hgccoord.lat, hgccoord.rad)
    return hgsframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, HelioCentric, HelioProjective)
def hcc_to_hpc(helioccoord, heliopframe):
    x = helioccoord.x.to(u.m)
    y = helioccoord.y.to(u.m)
    z = helioccoord.z.to(u.m)

    # d is calculated as the distance between the points
    # (x,y,z) and (0,0,D0).
    distance = np.sqrt(x**2 + y**2 + (helioccoord.D0.to(u.m) - z)**2)

    hpcx = np.rad2deg(np.arctan2(x, helioccoord.D0 - z))
    hpcy = np.rad2deg(np.arcsin(y / distance))

    representation = SphericalWrap180Representation(hpcx, hpcy, distance.to(u.km))
    return heliopframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, HelioProjective, HelioCentric)
def hpc_to_hcc(heliopcoord, heliocframe):
    heliopcoord = heliopcoord.calculate_distance()
    x = np.deg2rad(heliopcoord.Tx)
    y = np.deg2rad(heliopcoord.Ty)

    cosx = np.cos(x)
    sinx = np.sin(x)
    cosy = np.cos(y)
    siny = np.sin(y)

    rx = (heliopcoord.distance.to(u.m)) * cosy * sinx
    ry = (heliopcoord.distance.to(u.m)) * siny
    rz = (heliopcoord.D0.to(u.m)) - (heliopcoord.distance.to(u.m)) * cosy * cosx

    representation = CartesianRepresentation(rx.to(u.km), ry.to(u.km), rz.to(u.km))
    return heliocframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, HelioCentric, HelioGraphicStonyhurst)
def hcc_to_hgs(helioccoord, heliogframe):
    x = helioccoord.x.to(u.m)
    y = helioccoord.y.to(u.m)
    z = helioccoord.z.to(u.m)

    l0b0_pair = [helioccoord.L0, helioccoord.B0]

    l0_rad = l0b0_pair[0].to(u.rad)
    b0_deg = l0b0_pair[1]

    cosb = np.cos(np.deg2rad(b0_deg))
    sinb = np.sin(np.deg2rad(b0_deg))

    hecr = np.sqrt(x**2 + y**2 + z**2)
    hgln = np.arctan2(x, z * cosb - y * sinb) + l0_rad
    hglt = np.arcsin((y * cosb + z * sinb) / hecr)

    representation = SphericalWrap180Representation(np.rad2deg(hgln),
                                             np.rad2deg(hglt),
                                             hecr.to(u.km))
    return heliogframe.realize_frame(representation)


@frame_transform_graph.transform(FunctionTransform, HelioGraphicStonyhurst, HelioCentric)
def hgs_to_hcc(heliogcoord, heliocframe):
    hglon = heliogcoord.lon
    hglat = heliogcoord.lat
    r = heliogcoord.rad.to(u.m)

    l0b0_pair = [heliocframe.L0, heliocframe.B0]

    l0_rad = l0b0_pair[0].to(u.rad)
    b0_deg = l0b0_pair[1]

    lon = np.deg2rad(hglon)
    lat = np.deg2rad(hglat)

    cosb = np.cos(b0_deg.to(u.rad))
    sinb = np.sin(b0_deg.to(u.rad))

    lon = lon - l0_rad

    cosx = np.cos(lon)
    sinx = np.sin(lon)
    cosy = np.cos(lat)
    siny = np.sin(lat)

    x = r * cosy * sinx
    y = r * (siny * cosb - cosy * cosx * sinb)
    zz = r * (siny * sinb + cosy * cosx * cosb)

    representation = CartesianRepresentation(x.to(u.km), y.to(u.km), zz.to(u.km))
    return heliocframe.realize_frame(representation)
