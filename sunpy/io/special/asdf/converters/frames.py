from asdf_astropy.converters.coordinates.frame import FrameConverter


class SunpyFrameConverter(FrameConverter):
    def select_tag(self, obj, tags, ctx):
        # Sort the tags in reverse alphabetical order and pick the first (i.e.
        # the one with the highest version). This assumes that all the tags for
        # this converter are named the same other than the version number.
        tags = list(sorted(tags, reverse=True))
        return tags[0]


class HeliographicCarringtonConverter(SunpyFrameConverter):
    def from_yaml_tree(self, node, tag, ctx):
        # The 1.0.0 schema should be treated as having the observer at Earth
        if "1.0.0" in tag:
            node['frame_attributes']['observer'] = 'earth'
        return super().from_yaml_tree(node, tag, ctx)


SUNPY_FRAME_CONVERTERS = [
    SunpyFrameConverter(
        "tag:sunpy.org:sunpy/coordinates/frames/heliographic_stonyhurst-*",
        "sunpy.coordinates.frames.HeliographicStonyhurst"
    ),
    HeliographicCarringtonConverter(
        "tag:sunpy.org:sunpy/coordinates/frames/heliographic_carrington-*",
        "sunpy.coordinates.frames.HeliographicCarrington"
    ),
    SunpyFrameConverter(
        "tag:sunpy.org:sunpy/coordinates/frames/helioprojective-*",
        "sunpy.coordinates.frames.Helioprojective"
    ),
    SunpyFrameConverter(
        "tag:sunpy.org:sunpy/coordinates/frames/helioprojectiveradial-*",
        "sunpy.coordinates.frames.HelioprojectiveRadial"
    ),
    SunpyFrameConverter(
        "tag:sunpy.org:sunpy/coordinates/frames/heliocentricinertial-*",
        "sunpy.coordinates.frames.HeliocentricInertial"
    ),
    SunpyFrameConverter(
        "tag:sunpy.org:sunpy/coordinates/frames/heliocentricearthecliptic-*",
        "sunpy.coordinates.frames.HeliocentricEarthEcliptic"
    ),
    SunpyFrameConverter(
        "tag:sunpy.org:sunpy/coordinates/frames/heliocentric-*",
        "sunpy.coordinates.frames.Heliocentric"
    ),
    SunpyFrameConverter(
        "tag:sunpy.org:sunpy/coordinates/frames/geocentricsolarecliptic-*",
        "sunpy.coordinates.frames.GeocentricSolarEcliptic"
    ),
    SunpyFrameConverter(
        "tag:sunpy.org:sunpy/coordinates/frames/geocentricearthequatorial-*",
        "sunpy.coordinates.frames.GeocentricEarthEquatorial"
    ),
    SunpyFrameConverter(
        "tag:sunpy.org:sunpy/coordinates/frames/geomagnetic-*",
        "sunpy.coordinates.frames.Geomagnetic"
    ),
    SunpyFrameConverter(
        "tag:sunpy.org:sunpy/coordinates/frames/solarmagnetic-*",
        "sunpy.coordinates.frames.SolarMagnetic"
    ),
    SunpyFrameConverter(
        "tag:sunpy.org:sunpy/coordinates/frames/geocentricsolarmagnetospheric-*",
        "sunpy.coordinates.frames.GeocentricSolarMagnetospheric"
    ),
]
