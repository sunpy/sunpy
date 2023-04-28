# Author: Florian Mayer <florian.mayer@bitsrc.org>
#
# This module was developed with funding provided by
# the ESA Summer of Code (2011).

# The template can be found in tools/hektemplate.py

"""
This script is used to generate sunpy.net.hek.attrs. The rationale for using
code-generation in lieu of dynamic magic is that code-generation ensures
that tools (e.g. IPython) are able to automatically complete names of members.
Considering the sheer amount of members it is essential that users are able
to use completion.

Events are EventType objects. When they are directly ORed together, they are
joined together so that only one query is sent to the query. They may not
be ANDed together because an event cannot be of multiple types.

Events also have attributes which are _StringParamAttrWrapper, that means that
they overload the Python operators for strings and return a AttrComparison for
them. AttrComparison are used to specify the values of parameters of the event.
So, AR.NumSpots == 1 returns a AttrComparison that when encountered in a query
sets the GET parameters in a way that only active regions with only one spot
are returned. _StringParamAttrWrapper support <, <= ,>, >=, ==, != and like.
ComparisonParamAttrWrapper support all the operations mentioned above
barring like.
"""

# XXX: Maybe split into three modules and import them all into one so
# we do not need a template but generate one module in its entirety.

import os
import sys
from collections import defaultdict

EVENTS = [
    'AR', 'CME', 'CD', 'CH', 'CW', 'FI', 'FE', 'FA', 'FL', 'LP', 'OS', 'SS',
    'EF', 'CJ', 'PG', 'OT', 'NR', 'SG', 'SP', 'CR', 'CC', 'ER', 'TO'
]

# For some reason, the event type is "ce" but all its attributes start with
# "CME". This dict is here to consider this.
NAMES = defaultdict(lambda: None, {
    'CME': 'CE'
})
# These are just groups for attributes that are not _ListAttrs themselves.
OTHER = ['Area', 'BoundBox', 'Bound', 'OBS', 'Skel', 'FRM', 'Event', 'Outflow']
# There is no underscore after Wave in the names of the API, so we do not
# need to remove it.
OTHER_NOPAD = ['Wave', 'Veloc', 'Freq', 'Intens']
# Every attribute that neither starts with something in EVENTS, OTHER or
# OTHER_NOPAD, is put into the Misc class.

# XXX: Not all of them actually are string. We just use string for now because
# that is the type that has the most functionality.
fields = {
    'AR_CompactnessCls': '_StringParamAttrWrapper',
    'AR_IntensKurt': '_StringParamAttrWrapper',
    'AR_IntensMax': '_StringParamAttrWrapper',
    'AR_IntensMean': '_StringParamAttrWrapper',
    'AR_IntensMin': '_StringParamAttrWrapper',
    'AR_IntensSkew': '_StringParamAttrWrapper',
    'AR_IntensTotal': '_StringParamAttrWrapper',
    'AR_IntensUnit': '_StringParamAttrWrapper',
    'AR_IntensVar': '_StringParamAttrWrapper',
    'AR_McIntoshCls': '_StringParamAttrWrapper',
    'AR_MtWilsonCls': '_StringParamAttrWrapper',
    'AR_NOAANum': '_StringParamAttrWrapper',
    'AR_NOAAclass': '_StringParamAttrWrapper',
    'AR_NumSpots': '_StringParamAttrWrapper',
    'AR_PenumbraCls': '_StringParamAttrWrapper',
    'AR_Polarity': '_StringParamAttrWrapper',
    'AR_SpotAreaRaw': '_StringParamAttrWrapper',
    'AR_SpotAreaRawUncert': '_StringParamAttrWrapper',
    'AR_SpotAreaRawUnit': '_StringParamAttrWrapper',
    'AR_SpotAreaRepr': '_StringParamAttrWrapper',
    'AR_SpotAreaReprUncert': '_StringParamAttrWrapper',
    'AR_SpotAreaReprUnit': '_StringParamAttrWrapper',
    'AR_ZurichCls': '_StringParamAttrWrapper',
    'Area_AtDiskCenter': '_StringParamAttrWrapper',
    'Area_AtDiskCenterUncert': '_StringParamAttrWrapper',
    'Area_Raw': '_StringParamAttrWrapper',
    'Area_Uncert': '_StringParamAttrWrapper',
    'Area_Unit': '_StringParamAttrWrapper',
    'BoundBox_C1LL': '_StringParamAttrWrapper',
    'BoundBox_C1UR': '_StringParamAttrWrapper',
    'BoundBox_C2LL': '_StringParamAttrWrapper',
    'BoundBox_C2UR': '_StringParamAttrWrapper',
    'Bound_CCNsteps': '_StringParamAttrWrapper',
    'Bound_CCStartC1': '_StringParamAttrWrapper',
    'Bound_CCStartC2': '_StringParamAttrWrapper',
    'CC_AxisUnit': '_StringParamAttrWrapper',
    'CC_MajorAxis': '_StringParamAttrWrapper',
    'CC_MinorAxis': '_StringParamAttrWrapper',
    'CC_TiltAngleMajorFromRadial': '_StringParamAttrWrapper',
    'CC_TiltAngleUnit': '_StringParamAttrWrapper',
    'CD_Area': '_StringParamAttrWrapper',
    'CD_AreaUncert': '_StringParamAttrWrapper',
    'CD_AreaUnit': '_StringParamAttrWrapper',
    'CD_Mass': '_StringParamAttrWrapper',
    'CD_MassUncert': '_StringParamAttrWrapper',
    'CD_MassUnit': '_StringParamAttrWrapper',
    'CD_Volume': '_StringParamAttrWrapper',
    'CD_VolumeUncert': '_StringParamAttrWrapper',
    'CD_VolumeUnit': '_StringParamAttrWrapper',
    'CME_Accel': '_StringParamAttrWrapper',
    'CME_AccelUncert': '_StringParamAttrWrapper',
    'CME_AccelUnit': '_StringParamAttrWrapper',
    'CME_AngularWidth': '_StringParamAttrWrapper',
    'CME_AngularWidthUnit': '_StringParamAttrWrapper',
    'CME_Mass': '_StringParamAttrWrapper',
    'CME_MassUncert': '_StringParamAttrWrapper',
    'CME_MassUnit': '_StringParamAttrWrapper',
    'CME_RadialLinVel': '_StringParamAttrWrapper',
    'CME_RadialLinVelMax': '_StringParamAttrWrapper',
    'CME_RadialLinVelMin': '_StringParamAttrWrapper',
    'CME_RadialLinVelStddev': '_StringParamAttrWrapper',
    'CME_RadialLinVelUncert': '_StringParamAttrWrapper',
    'CME_RadialLinVelUnit': '_StringParamAttrWrapper',
    'EF_AspectRatio': '_StringParamAttrWrapper',
    'EF_AxisLength': '_StringParamAttrWrapper',
    'EF_AxisOrientation': '_StringParamAttrWrapper',
    'EF_AxisOrientationUnit': '_StringParamAttrWrapper',
    'EF_FluxUnit': '_StringParamAttrWrapper',
    'EF_LengthUnit': '_StringParamAttrWrapper',
    'EF_NegEquivRadius': '_StringParamAttrWrapper',
    'EF_NegPeakFluxOnsetRate': '_StringParamAttrWrapper',
    'EF_OnsetRateUnit': '_StringParamAttrWrapper',
    'EF_PosEquivRadius': '_StringParamAttrWrapper',
    'EF_PosPeakFluxOnsetRate': '_StringParamAttrWrapper',
    'EF_ProximityRatio': '_StringParamAttrWrapper',
    'EF_SumNegSignedFlux': '_StringParamAttrWrapper',
    'EF_SumPosSignedFlux': '_StringParamAttrWrapper',
    'Event_C1Error': '_StringParamAttrWrapper',
    'Event_C2Error': '_StringParamAttrWrapper',
    'Event_ClippedSpatial': '_StringParamAttrWrapper',
    'Event_ClippedTemporal': '_StringParamAttrWrapper',
    'Event_Coord1': '_StringParamAttrWrapper',
    'Event_Coord2': '_StringParamAttrWrapper',
    'Event_Coord3': '_StringParamAttrWrapper',
    'Event_CoordSys': '_StringParamAttrWrapper',
    'Event_CoordUnit': '_StringParamAttrWrapper',
    'Event_MapURL': '_StringParamAttrWrapper',
    'Event_MaskURL': '_StringParamAttrWrapper',
    'Event_Npixels': '_StringParamAttrWrapper',
    'Event_PixelUnit': '_StringParamAttrWrapper',
    'Event_Probability': '_StringParamAttrWrapper',
    'Event_TestFlag': '_StringParamAttrWrapper',
    'Event_Type': '_StringParamAttrWrapper',
    'FI_BarbsL': '_StringParamAttrWrapper',
    'FI_BarbsR': '_StringParamAttrWrapper',
    'FI_BarbsTot': '_StringParamAttrWrapper',
    'FI_Chirality': '_StringParamAttrWrapper',
    'FI_Length': '_StringParamAttrWrapper',
    'FI_LengthUnit': '_StringParamAttrWrapper',
    'FI_Tilt': '_StringParamAttrWrapper',
    'FL_EFoldTime': '_StringParamAttrWrapper',
    'FL_EFoldTimeUnit': '_StringParamAttrWrapper',
    'FL_Fluence': '_StringParamAttrWrapper',
    'FL_FluenceUnit': '_StringParamAttrWrapper',
    'FL_GOESCls': '_StringParamAttrWrapper',
    'FL_PeakEM': '_StringParamAttrWrapper',
    'FL_PeakEMUnit': '_StringParamAttrWrapper',
    'FL_PeakFlux': '_StringParamAttrWrapper',
    'FL_PeakFluxUnit': '_StringParamAttrWrapper',
    'FL_PeakTemp': '_StringParamAttrWrapper',
    'FL_PeakTempUnit': '_StringParamAttrWrapper',
    'FRM_Contact': '_StringParamAttrWrapper',
    'FRM_HumanFlag': '_StringParamAttrWrapper',
    'FRM_Identifier': '_StringParamAttrWrapper',
    'FRM_Institute': '_StringParamAttrWrapper',
    'FRM_Name': '_StringParamAttrWrapper',
    'FRM_ParamSet': '_StringParamAttrWrapper',
    'FRM_SpecificID': '_StringParamAttrWrapper',
    'FRM_URL': '_StringParamAttrWrapper',
    'FRM_VersionNumber': '_StringParamAttrWrapper',
    'FreqMaxRange': '_StringParamAttrWrapper',
    'FreqMinRange': '_StringParamAttrWrapper',
    'FreqPeakPower': '_StringParamAttrWrapper',
    'FreqUnit': '_StringParamAttrWrapper',
    'IntensMaxAmpl': '_StringParamAttrWrapper',
    'IntensMinAmpl': '_StringParamAttrWrapper',
    'IntensUnit': '_StringParamAttrWrapper',
    'KB_Archivist': '_StringParamAttrWrapper',
    'MaxMagFieldStrength': '_StringParamAttrWrapper',
    'MaxMagFieldStrengthUnit': '_StringParamAttrWrapper',
    'OBS_ChannelID': '_StringParamAttrWrapper',
    'OBS_DataPrepURL': '_StringParamAttrWrapper',
    'OBS_FirstProcessingDate': '_StringParamAttrWrapper',
    'OBS_IncludesNRT': '_StringParamAttrWrapper',
    'OBS_Instrument': '_StringParamAttrWrapper',
    'OBS_LastProcessingDate': '_StringParamAttrWrapper',
    'OBS_LevelNum': '_StringParamAttrWrapper',
    'OBS_MeanWavel': '_StringParamAttrWrapper',
    'OBS_Observatory': '_StringParamAttrWrapper',
    'OBS_Title': '_StringParamAttrWrapper',
    'OBS_WavelUnit': '_StringParamAttrWrapper',
    'OscillNPeriods': '_StringParamAttrWrapper',
    'OscillNPeriodsUncert': '_StringParamAttrWrapper',
    'Outflow_Length': '_StringParamAttrWrapper',
    'Outflow_LengthUnit': '_StringParamAttrWrapper',
    'Outflow_OpeningAngle': '_StringParamAttrWrapper',
    'Outflow_Speed': '_StringParamAttrWrapper',
    'Outflow_SpeedUnit': '_StringParamAttrWrapper',
    'Outflow_TransSpeed': '_StringParamAttrWrapper',
    'Outflow_Width': '_StringParamAttrWrapper',
    'Outflow_WidthUnit': '_StringParamAttrWrapper',
    'PeakPower': '_StringParamAttrWrapper',
    'PeakPowerUnit': '_StringParamAttrWrapper',
    'RasterScanType': '_StringParamAttrWrapper',
    'SG_AspectRatio': '_StringParamAttrWrapper',
    'SG_Chirality': '_StringParamAttrWrapper',
    'SG_MeanContrast': '_StringParamAttrWrapper',
    'SG_Orientation': '_StringParamAttrWrapper',
    'SG_PeakContrast': '_StringParamAttrWrapper',
    'SG_Shape': '_StringParamAttrWrapper',
    'SS_SpinRate': '_StringParamAttrWrapper',
    'SS_SpinRateUnit': '_StringParamAttrWrapper',
    'Skel_Curvature': '_StringParamAttrWrapper',
    'Skel_Nsteps': '_StringParamAttrWrapper',
    'Skel_StartC1': '_StringParamAttrWrapper',
    'Skel_StartC2': '_StringParamAttrWrapper',
    'TO_Shape': '_StringParamAttrWrapper',
    'VelocMaxAmpl': '_StringParamAttrWrapper',
    'VelocMaxPower': '_StringParamAttrWrapper',
    'VelocMaxPowerUncert': '_StringParamAttrWrapper',
    'VelocMinAmpl': '_StringParamAttrWrapper',
    'VelocUnit': '_StringParamAttrWrapper',
    'WaveDisplMaxAmpl': '_StringParamAttrWrapper',
    'WaveDisplMinAmpl': '_StringParamAttrWrapper',
    'WaveDisplUnit': '_StringParamAttrWrapper',
    'WavelMaxPower': '_StringParamAttrWrapper',
    'WavelMaxPowerUncert': '_StringParamAttrWrapper',
    'WavelMaxRange': '_StringParamAttrWrapper',
    'WavelMinRange': '_StringParamAttrWrapper',
    'WavelUnit': '_StringParamAttrWrapper'
}


def mk_gen(rest):
    """ Generate Misc class. """
    ret = ''
    ret += '@_makeinstance\nclass Misc:\n'
    for elem in sorted(rest):
        ret += f'    {elem} = {fields[elem]}({elem!r})\n'
    return ret


def mk_cls(key, used, pad=1, nokeys=True, init=True, name=None, base='EventType'):
    if name is None:
        name = key

    keys = sorted(
        [(k, v) for k, v in fields.items() if k.startswith(key)]
    )
    used.update({k for k, v in keys})
    if not keys:
        if not nokeys:
            raise ValueError
        return f'{key} = EventType({name.lower()!r})'
    ret = ''
    if base != 'object':
        ret += f'@_makeinstance\nclass {name}({base}):\n'
    else:
        ret += '@_makeinstance\nclass %s:\n' % name
    for k, v in keys:
        ret += f'    {k[len(key) + pad:]} = {v}({k!r})\n'
    if init:
        ret += '''\n    def __init__(self):
        super().__init__(%r)''' % name.lower()
    return ret


if __name__ == '__main__':
    BUFFER = 4096
    used = set()
    tmpl = (
        os.path.join(os.path.dirname(__file__), 'hektemplate.py')
        if len(sys.argv) <= 2 else sys.argv[2]
    )
    dest = (
        os.path.join(
            os.path.dirname(__file__), os.pardir, 'sunpy', 'net', 'hek',
            'attrs.py')
        if len(sys.argv) <= 1 else sys.argv[1]
    )

    if dest == '-':
        fd = sys.stdout
    else:
        fd = open(dest, 'w')

    tmplfd = open(tmpl)

    while True:
        buf = tmplfd.read(BUFFER)
        if not buf:
            break
        fd.write(buf)

    fd.write('\n\n')
    fd.write('\n\n\n'.join(mk_cls(evt, used, name=NAMES[evt]) for evt in EVENTS))
    fd.write('\n\n\n')
    fd.write('\n\n'.join(mk_cls(evt, used, 0, 0, 0, NAMES[evt], 'object') for evt in OTHER_NOPAD))
    fd.write('\n\n')
    fd.write('\n\n'.join(mk_cls(evt, used, 1, 0, 0, NAMES[evt], 'object') for evt in OTHER))
    fd.write('\n\n')
    fd.write(mk_gen(set(fields) - used))
