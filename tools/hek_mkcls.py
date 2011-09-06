import shutil
import sys
import os

from collections import defaultdict

EVENTS = [
    'AR', 'CME', 'CD', 'CH', 'CW', 'FI', 'FE', 'FA', 'FL', 'LP', 'OS', 'SS',
    'EF', 'CJ', 'PG', 'OT', 'NR', 'SG', 'SP', 'CR', 'CC', 'ER', 'TO'
]

NAMES = defaultdict(lambda: None, {
    'CME': 'CE'
})

OTHER = ['Area', 'BoundBox', 'Bound', 'OBS', 'Skel', 'FRM', 'Event', 'Outflow']
OTHER_NOPAD = ['Wave', 'Veloc', 'Freq', 'Intens']

fields = {
    'AR_CompactnessCls': 'StringParamAttrWrapper',
    'AR_IntensKurt': 'StringParamAttrWrapper',
    'AR_IntensMax': 'StringParamAttrWrapper',
    'AR_IntensMean': 'StringParamAttrWrapper',
    'AR_IntensMin': 'StringParamAttrWrapper',
    'AR_IntensSkew': 'StringParamAttrWrapper',
    'AR_IntensTotal': 'StringParamAttrWrapper',
    'AR_IntensUnit': 'StringParamAttrWrapper',
    'AR_IntensVar': 'StringParamAttrWrapper',
    'AR_McIntoshCls': 'StringParamAttrWrapper',
    'AR_MtWilsonCls': 'StringParamAttrWrapper',
    'AR_NOAANum': 'StringParamAttrWrapper',
    'AR_NOAAclass': 'StringParamAttrWrapper',
    'AR_NumSpots': 'StringParamAttrWrapper',
    'AR_PenumbraCls': 'StringParamAttrWrapper',
    'AR_Polarity': 'StringParamAttrWrapper',
    'AR_SpotAreaRaw': 'StringParamAttrWrapper',
    'AR_SpotAreaRawUncert': 'StringParamAttrWrapper',
    'AR_SpotAreaRawUnit': 'StringParamAttrWrapper',
    'AR_SpotAreaRepr': 'StringParamAttrWrapper',
    'AR_SpotAreaReprUncert': 'StringParamAttrWrapper',
    'AR_SpotAreaReprUnit': 'StringParamAttrWrapper',
    'AR_ZurichCls': 'StringParamAttrWrapper',
    'Area_AtDiskCenter': 'StringParamAttrWrapper',
    'Area_AtDiskCenterUncert': 'StringParamAttrWrapper',
    'Area_Raw': 'StringParamAttrWrapper',
    'Area_Uncert': 'StringParamAttrWrapper',
    'Area_Unit': 'StringParamAttrWrapper',
    'BoundBox_C1LL': 'StringParamAttrWrapper',
    'BoundBox_C1UR': 'StringParamAttrWrapper',
    'BoundBox_C2LL': 'StringParamAttrWrapper',
    'BoundBox_C2UR': 'StringParamAttrWrapper',
    'Bound_CCNsteps': 'StringParamAttrWrapper',
    'Bound_CCStartC1': 'StringParamAttrWrapper',
    'Bound_CCStartC2': 'StringParamAttrWrapper',
    'CC_AxisUnit': 'StringParamAttrWrapper',
    'CC_MajorAxis': 'StringParamAttrWrapper',
    'CC_MinorAxis': 'StringParamAttrWrapper',
    'CC_TiltAngleMajorFromRadial': 'StringParamAttrWrapper',
    'CC_TiltAngleUnit': 'StringParamAttrWrapper',
    'CD_Area': 'StringParamAttrWrapper',
    'CD_AreaUncert': 'StringParamAttrWrapper',
    'CD_AreaUnit': 'StringParamAttrWrapper',
    'CD_Mass': 'StringParamAttrWrapper',
    'CD_MassUncert': 'StringParamAttrWrapper',
    'CD_MassUnit': 'StringParamAttrWrapper',
    'CD_Volume': 'StringParamAttrWrapper',
    'CD_VolumeUncert': 'StringParamAttrWrapper',
    'CD_VolumeUnit': 'StringParamAttrWrapper',
    'CME_Accel': 'StringParamAttrWrapper',
    'CME_AccelUncert': 'StringParamAttrWrapper',
    'CME_AccelUnit': 'StringParamAttrWrapper',
    'CME_AngularWidth': 'StringParamAttrWrapper',
    'CME_AngularWidthUnit': 'StringParamAttrWrapper',
    'CME_Mass': 'StringParamAttrWrapper',
    'CME_MassUncert': 'StringParamAttrWrapper',
    'CME_MassUnit': 'StringParamAttrWrapper',
    'CME_RadialLinVel': 'StringParamAttrWrapper',
    'CME_RadialLinVelMax': 'StringParamAttrWrapper',
    'CME_RadialLinVelMin': 'StringParamAttrWrapper',
    'CME_RadialLinVelStddev': 'StringParamAttrWrapper',
    'CME_RadialLinVelUncert': 'StringParamAttrWrapper',
    'CME_RadialLinVelUnit': 'StringParamAttrWrapper',
    'EF_AspectRatio': 'StringParamAttrWrapper',
    'EF_AxisLength': 'StringParamAttrWrapper',
    'EF_AxisOrientation': 'StringParamAttrWrapper',
    'EF_AxisOrientationUnit': 'StringParamAttrWrapper',
    'EF_FluxUnit': 'StringParamAttrWrapper',
    'EF_LengthUnit': 'StringParamAttrWrapper',
    'EF_NegEquivRadius': 'StringParamAttrWrapper',
    'EF_NegPeakFluxOnsetRate': 'StringParamAttrWrapper',
    'EF_OnsetRateUnit': 'StringParamAttrWrapper',
    'EF_PosEquivRadius': 'StringParamAttrWrapper',
    'EF_PosPeakFluxOnsetRate': 'StringParamAttrWrapper',
    'EF_ProximityRatio': 'StringParamAttrWrapper',
    'EF_SumNegSignedFlux': 'StringParamAttrWrapper',
    'EF_SumPosSignedFlux': 'StringParamAttrWrapper',
    'Event_C1Error': 'StringParamAttrWrapper',
    'Event_C2Error': 'StringParamAttrWrapper',
    'Event_ClippedSpatial': 'StringParamAttrWrapper',
    'Event_ClippedTemporal': 'StringParamAttrWrapper',
    'Event_Coord1': 'StringParamAttrWrapper',
    'Event_Coord2': 'StringParamAttrWrapper',
    'Event_Coord3': 'StringParamAttrWrapper',
    'Event_CoordSys': 'StringParamAttrWrapper',
    'Event_CoordUnit': 'StringParamAttrWrapper',
    'Event_MapURL': 'StringParamAttrWrapper',
    'Event_MaskURL': 'StringParamAttrWrapper',
    'Event_Npixels': 'StringParamAttrWrapper',
    'Event_PixelUnit': 'StringParamAttrWrapper',
    'Event_Probability': 'StringParamAttrWrapper',
    'Event_TestFlag': 'StringParamAttrWrapper',
    'Event_Type': 'StringParamAttrWrapper',
    'FI_BarbsL': 'StringParamAttrWrapper',
    'FI_BarbsR': 'StringParamAttrWrapper',
    'FI_BarbsTot': 'StringParamAttrWrapper',
    'FI_Chirality': 'StringParamAttrWrapper',
    'FI_Length': 'StringParamAttrWrapper',
    'FI_LengthUnit': 'StringParamAttrWrapper',
    'FI_Tilt': 'StringParamAttrWrapper',
    'FL_EFoldTime': 'StringParamAttrWrapper',
    'FL_EFoldTimeUnit': 'StringParamAttrWrapper',
    'FL_Fluence': 'StringParamAttrWrapper',
    'FL_FluenceUnit': 'StringParamAttrWrapper',
    'FL_GOESCls': 'StringParamAttrWrapper',
    'FL_PeakEM': 'StringParamAttrWrapper',
    'FL_PeakEMUnit': 'StringParamAttrWrapper',
    'FL_PeakFlux': 'StringParamAttrWrapper',
    'FL_PeakFluxUnit': 'StringParamAttrWrapper',
    'FL_PeakTemp': 'StringParamAttrWrapper',
    'FL_PeakTempUnit': 'StringParamAttrWrapper',
    'FRM_Contact': 'StringParamAttrWrapper',
    'FRM_HumanFlag': 'StringParamAttrWrapper',
    'FRM_Identifier': 'StringParamAttrWrapper',
    'FRM_Institute': 'StringParamAttrWrapper',
    'FRM_Name': 'StringParamAttrWrapper',
    'FRM_ParamSet': 'StringParamAttrWrapper',
    'FRM_SpecificID': 'StringParamAttrWrapper',
    'FRM_URL': 'StringParamAttrWrapper',
    'FRM_VersionNumber': 'StringParamAttrWrapper',
    'FreqMaxRange': 'StringParamAttrWrapper',
    'FreqMinRange': 'StringParamAttrWrapper',
    'FreqPeakPower': 'StringParamAttrWrapper',
    'FreqUnit': 'StringParamAttrWrapper',
    'IntensMaxAmpl': 'StringParamAttrWrapper',
    'IntensMinAmpl': 'StringParamAttrWrapper',
    'IntensUnit': 'StringParamAttrWrapper',
    'KB_Archivist': 'StringParamAttrWrapper',
    'MaxMagFieldStrength': 'StringParamAttrWrapper',
    'MaxMagFieldStrengthUnit': 'StringParamAttrWrapper',
    'OBS_ChannelID': 'StringParamAttrWrapper',
    'OBS_DataPrepURL': 'StringParamAttrWrapper',
    'OBS_FirstProcessingDate': 'StringParamAttrWrapper',
    'OBS_IncludesNRT': 'StringParamAttrWrapper',
    'OBS_Instrument': 'StringParamAttrWrapper',
    'OBS_LastProcessingDate': 'StringParamAttrWrapper',
    'OBS_LevelNum': 'StringParamAttrWrapper',
    'OBS_MeanWavel': 'StringParamAttrWrapper',
    'OBS_Observatory': 'StringParamAttrWrapper',
    'OBS_Title': 'StringParamAttrWrapper',
    'OBS_WavelUnit': 'StringParamAttrWrapper',
    'OscillNPeriods': 'StringParamAttrWrapper',
    'OscillNPeriodsUncert': 'StringParamAttrWrapper',
    'Outflow_Length': 'StringParamAttrWrapper',
    'Outflow_LengthUnit': 'StringParamAttrWrapper',
    'Outflow_OpeningAngle': 'StringParamAttrWrapper',
    'Outflow_Speed': 'StringParamAttrWrapper',
    'Outflow_SpeedUnit': 'StringParamAttrWrapper',
    'Outflow_TransSpeed': 'StringParamAttrWrapper',
    'Outflow_Width': 'StringParamAttrWrapper',
    'Outflow_WidthUnit': 'StringParamAttrWrapper',
    'PeakPower': 'StringParamAttrWrapper',
    'PeakPowerUnit': 'StringParamAttrWrapper',
    'RasterScanType': 'StringParamAttrWrapper',
    'SG_AspectRatio': 'StringParamAttrWrapper',
    'SG_Chirality': 'StringParamAttrWrapper',
    'SG_MeanContrast': 'StringParamAttrWrapper',
    'SG_Orientation': 'StringParamAttrWrapper',
    'SG_PeakContrast': 'StringParamAttrWrapper',
    'SG_Shape': 'StringParamAttrWrapper',
    'SS_SpinRate': 'StringParamAttrWrapper',
    'SS_SpinRateUnit': 'StringParamAttrWrapper',
    'Skel_Curvature': 'StringParamAttrWrapper',
    'Skel_Nsteps': 'StringParamAttrWrapper',
    'Skel_StartC1': 'StringParamAttrWrapper',
    'Skel_StartC2': 'StringParamAttrWrapper',
    'TO_Shape': 'StringParamAttrWrapper',
    'VelocMaxAmpl': 'StringParamAttrWrapper',
    'VelocMaxPower': 'StringParamAttrWrapper',
    'VelocMaxPowerUncert': 'StringParamAttrWrapper',
    'VelocMinAmpl': 'StringParamAttrWrapper',
    'VelocUnit': 'StringParamAttrWrapper',
    'WaveDisplMaxAmpl': 'StringParamAttrWrapper',
    'WaveDisplMinAmpl': 'StringParamAttrWrapper',
    'WaveDisplUnit': 'StringParamAttrWrapper',
    'WavelMaxPower': 'StringParamAttrWrapper',
    'WavelMaxPowerUncert': 'StringParamAttrWrapper',
    'WavelMaxRange': 'StringParamAttrWrapper',
    'WavelMinRange': 'StringParamAttrWrapper',
    'WavelUnit': 'StringParamAttrWrapper'
}

def mk_gen(rest):
    ret = ''
    ret += 'class Misc(object):\n'
    for elem in sorted(rest):
        ret += '    %s = %s(%r)\n' %(elem, fields[elem], elem)
    return ret

def mk_cls(key, used, pad=1, nokeys=True, init=True, name=None):
    if name is None:
        name = key
    
    keys = sorted(
        [(k, v) for k, v in fields.iteritems() if k.startswith(key)]
    )
    used.update(set([k for k, v in keys]))
    if not keys:
        if not nokeys:
            raise ValueError
        return '%s = ListAttr("event_type", %r)' % (key, name.lower())
    ret = ''
    ret += '@apply\nclass %s(ListAttr):\n' % name
    for k, v in keys:
        ret += '    %s = %s(%r)\n' % (k[len(key) + pad:], v, k)
    if init:
        ret += '''    def __init__(self):
            ListAttr.__init__(self, "event_type", %r)''' % name.lower()
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
    
    fd.write(
        "# THIS IS AN AUTOMATICALLY GENERATED FILE\n"
        "# DO NOT EDIT!\n"
    )
    
    while True:
        buf = tmplfd.read(BUFFER)
        if not buf:
            break
        fd.write(buf)
    
    fd.write('\n\n')
    fd.write('\n\n'.join(mk_cls(evt, used, name=NAMES[evt]) for evt in EVENTS))
    fd.write('\n\n')
    fd.write('\n\n'.join(mk_cls(evt, used, 0, 0, 0, name=NAMES[evt]) for evt in OTHER_NOPAD))
    fd.write('\n\n')
    fd.write('\n\n'.join(mk_cls(evt, used, 1, 0, 0, name=NAMES[evt]) for evt in OTHER))
    fd.write('\n\n')
    fd.write(mk_gen(set(fields) - used))
    

