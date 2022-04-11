from asdf_astropy.converters.coordinates.frame import FrameConverter


# HeliographicCarrington has multiple schema versions
class HeliographicCarringtonConverter(FrameConverter):
    def from_yaml_tree(self, node, tag, ctx):
        # The 1.0.0 schema should be treated as having the observer at Earth
        if "1.0.0" in tag:
            node['frame_attributes']['observer'] = 'earth'
        return super().from_yaml_tree(node, tag, ctx)


SUNPY_FRAME_CONVERTERS = [
    FrameConverter(
        "tag:sunpy.org/sunpy/frames/heliographic_stonyhurst-*",
        "sunpy.coordinates.frames.HeliographicStonyhurst"
    ),
    HeliographicCarringtonConverter(
        "tag:sunpy.org/sunpy/frames/heliographic_carrington-*",
        "sunpy.coordinates.frames.HeliographicCarrington"
    ),
    FrameConverter(
        "tag:sunpy.org/sunpy/frames/helioprojective-*",
        "sunpy.coordinates.frames.Helioprojective"
    ),
    FrameConverter(
        "tag:sunpy.org/sunpy/frames/heliocentricinertial-*",
        "sunpy.coordinates.frames.HeliocentricInertial"
    ),
    FrameConverter(
        "tag:sunpy.org/sunpy/frames/heliocentricearthecliptic-*",
        "sunpy.coordinates.frames.HeliocentricEarthEcliptic"
    ),
    FrameConverter(
        "tag:sunpy.org/sunpy/frames/heliocentric-*",
        "sunpy.coordinates.frames.Heliocentric"
    ),
    FrameConverter(
        "tag:sunpy.org/sunpy/frames/geocentricsolarecliptic-*",
        "sunpy.coordinates.frames.GeocentricSolarEcliptic"
    ),
    FrameConverter(
        "tag:sunpy.org/sunpy/frames/geocentricearthequatorial-*",
        "sunpy.coordinates.frames.GeocentricEarthEquatorial"
    ),
]
