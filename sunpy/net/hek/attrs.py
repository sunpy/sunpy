# THIS IS AN AUTOMATICALLY GENERATED FILE
# DO NOT EDIT!
# The template can be found in tools/hektemplate.py

from datetime import datetime
from sunpy.net import attr
from sunpy.util import util

class ParamAttr(attr.Attr):
    def __init__(self, name, op, value):
        attr.Attr.__init__(self)
        self.name = name
        self.op = op
        self.value = value
    
    def collides(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.op == other.op and self.name == other.name


class BoolParamAttr(ParamAttr):
    def __init__(self, name, value='true'):
        ParamAttr.__init__(self, name, '=', value)
    
    def __neg__(self):
        if self.value == 'true':
            return BoolParamAttr(self.name, 'false')
        else:
            return BoolParamAttr(self.name)
    
    def __pos__(self):
        return BoolParamAttr(self.name)


class ListAttr(attr.Attr):
    def __init__(self, key, item):
        attr.Attr.__init__(self)
        
        self.key = key
        self.item = item
    
    def collides(self, other):
        return False
    
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return vars(self) == vars(other)
    
    def __hash__(self):
        return hash(tuple(vars(self).itervalues()))


class EventType(ListAttr):
    def __init__(self, item):
        ListAttr.__init__(self, 'event_type', item)


class Time(attr.Attr):
    def __init__(self, start, end):
        attr.Attr.__init__(self)
        self.start = start
        self.end = end
    
    def collides(self, other):
        return isinstance(other, Time)
    
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return vars(self) == vars(other)
    
    def __hash__(self):
        return hash(tuple(vars(self).itervalues()))
    
    @classmethod
    def dt(cls, start, end):
        return cls(datetime(*start), datetime(*end))


# pylint: disable=R0913
class SpartialRegion(attr.Attr):
    def __init__(
        self, x1=-1200, y1=-1200, x2=1200, y2=1200, sys='helioprojective'):
        attr.Attr.__init__(self)
        
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.sys = sys
    
    def collides(self, other):
        return isinstance(other, SpartialRegion)
    
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return vars(self) == vars(other)
    
    def __hash__(self):
        return hash(tuple(vars(self).itervalues()))


class Contains(attr.Attr):
    def __init__(self, *types):
        attr.Attr.__init__(self)
        self.types = types
    
    def collides(self, other):
        return False
    
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return vars(self) == vars(other)
    
    def __hash__(self):
        return hash(tuple(vars(self).itervalues()))


class ComparisonParamAttrWrapper(object):
    def __init__(self, name):
        self.name = name
    
    def __lt__(self, other):
        return ParamAttr(self.name, '<', other)
    
    def __le__(self, other):
        return ParamAttr(self.name, '<=', other)
    
    def __gt__(self, other):
        return ParamAttr(self.name, '>', other)
    
    def __ge__(self, other):
        return ParamAttr(self.name, '>=', other)
    
    def __eq__(self, other):
        return ParamAttr(self.name, '=', other)
    
    def __neq__(self, other):
        return ParamAttr(self.name, '!=', other)


class StringParamAttrWrapper(ComparisonParamAttrWrapper):
    def like(self, other):
        return ParamAttr(self.name, 'like', other)


class NumberParamAttrWrapper(ComparisonParamAttrWrapper):
    pass


walker = attr.AttrWalker()

@walker.add_applier(Contains)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    dct['type'] = 'contains'
    if not Contains in state:
        state[Contains] = 1
    
    nid = state[Contains]
    n = 0
    for n, type_ in enumerate(root.types):
        dct['event_type%d' % (nid + n)] = type_
    state[Contains] += n
    return dct

@walker.add_creator(
    Time, SpartialRegion, ListAttr, ParamAttr, attr.AttrAnd, Contains)
# pylint: disable=E0102,C0103,W0613
def _c(wlk, root, state):
    value = {}
    wlk.apply(root, state, value)
    return [value]

@walker.add_applier(Time)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    dct['event_starttime'] = util.anytim(root.start).strftime('%Y-%m-%dT%H:%M:%S')
    dct['event_endtime'] = util.anytim(root.end).strftime('%Y-%m-%dT%H:%M:%S')
    return dct

@walker.add_applier(SpartialRegion)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    dct['x1'] = root.x1
    dct['y1'] = root.y1
    dct['x2'] = root.x2
    dct['y2'] = root.y2
    dct['event_coordsys'] = root.sys
    return dct

@walker.add_applier(ListAttr)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    if root.key in dct:
        dct[root.key] += ',%s' % root.item
    else:
        dct[root.key] = root.item
    return dct

@walker.add_applier(EventType)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    if dct.get('type', None) == 'contains':
        raise ValueError
    
    return wlk.super_apply(super(EventType, root), state, dct)

@walker.add_applier(ParamAttr)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    if not ParamAttr in state:
        state[ParamAttr] = 0
    
    nid = state[ParamAttr]
    dct['param%d' % nid] = root.name
    dct['op%d' % nid] = root.op
    dct['value%d' % nid] = root.value
    state[ParamAttr] += 1
    return dct

@walker.add_applier(attr.AttrAnd)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    for attribute in root.attrs:
        wlk.apply(attribute, state, dct)

@walker.add_creator(attr.AttrOr)
# pylint: disable=E0102,C0103,W0613
def _c(wlk, root, state):
    blocks = []
    for attribute in root.attrs:
        blocks.extend(wlk.create(attribute, state))
    return blocks

@walker.add_creator(attr.DummyAttr)
# pylint: disable=E0102,C0103,W0613
def _c(wlk, root, state):
    return {}

@walker.add_applier(attr.DummyAttr)
# pylint: disable=E0102,C0103,W0613
def _a(wlk, root, state, dct):
    pass


@apply
class AR(ListAttr):
    CompactnessCls = StringParamAttrWrapper('AR_CompactnessCls')
    IntensKurt = StringParamAttrWrapper('AR_IntensKurt')
    IntensMax = StringParamAttrWrapper('AR_IntensMax')
    IntensMean = StringParamAttrWrapper('AR_IntensMean')
    IntensMin = StringParamAttrWrapper('AR_IntensMin')
    IntensSkew = StringParamAttrWrapper('AR_IntensSkew')
    IntensTotal = StringParamAttrWrapper('AR_IntensTotal')
    IntensUnit = StringParamAttrWrapper('AR_IntensUnit')
    IntensVar = StringParamAttrWrapper('AR_IntensVar')
    McIntoshCls = StringParamAttrWrapper('AR_McIntoshCls')
    MtWilsonCls = StringParamAttrWrapper('AR_MtWilsonCls')
    NOAANum = StringParamAttrWrapper('AR_NOAANum')
    NOAAclass = StringParamAttrWrapper('AR_NOAAclass')
    NumSpots = StringParamAttrWrapper('AR_NumSpots')
    PenumbraCls = StringParamAttrWrapper('AR_PenumbraCls')
    Polarity = StringParamAttrWrapper('AR_Polarity')
    SpotAreaRaw = StringParamAttrWrapper('AR_SpotAreaRaw')
    SpotAreaRawUncert = StringParamAttrWrapper('AR_SpotAreaRawUncert')
    SpotAreaRawUnit = StringParamAttrWrapper('AR_SpotAreaRawUnit')
    SpotAreaRepr = StringParamAttrWrapper('AR_SpotAreaRepr')
    SpotAreaReprUncert = StringParamAttrWrapper('AR_SpotAreaReprUncert')
    SpotAreaReprUnit = StringParamAttrWrapper('AR_SpotAreaReprUnit')
    ZurichCls = StringParamAttrWrapper('AR_ZurichCls')
    def __init__(self):
            ListAttr.__init__(self, "event_type", 'ar')

CE = ListAttr("event_type", 'ce')

@apply
class CD(ListAttr):
    Area = StringParamAttrWrapper('CD_Area')
    AreaUncert = StringParamAttrWrapper('CD_AreaUncert')
    AreaUnit = StringParamAttrWrapper('CD_AreaUnit')
    Mass = StringParamAttrWrapper('CD_Mass')
    MassUncert = StringParamAttrWrapper('CD_MassUncert')
    MassUnit = StringParamAttrWrapper('CD_MassUnit')
    Volume = StringParamAttrWrapper('CD_Volume')
    VolumeUncert = StringParamAttrWrapper('CD_VolumeUncert')
    VolumeUnit = StringParamAttrWrapper('CD_VolumeUnit')
    def __init__(self):
            ListAttr.__init__(self, "event_type", 'cd')

CH = ListAttr("event_type", 'ch')

CW = ListAttr("event_type", 'cw')

@apply
class FI(ListAttr):
    BarbsL = StringParamAttrWrapper('FI_BarbsL')
    BarbsR = StringParamAttrWrapper('FI_BarbsR')
    BarbsTot = StringParamAttrWrapper('FI_BarbsTot')
    Chirality = StringParamAttrWrapper('FI_Chirality')
    Length = StringParamAttrWrapper('FI_Length')
    LengthUnit = StringParamAttrWrapper('FI_LengthUnit')
    Tilt = StringParamAttrWrapper('FI_Tilt')
    def __init__(self):
            ListAttr.__init__(self, "event_type", 'fi')

FE = ListAttr("event_type", 'fe')

FA = ListAttr("event_type", 'fa')

@apply
class FL(ListAttr):
    EFoldTime = StringParamAttrWrapper('FL_EFoldTime')
    EFoldTimeUnit = StringParamAttrWrapper('FL_EFoldTimeUnit')
    Fluence = StringParamAttrWrapper('FL_Fluence')
    FluenceUnit = StringParamAttrWrapper('FL_FluenceUnit')
    GOESCls = StringParamAttrWrapper('FL_GOESCls')
    PeakEM = StringParamAttrWrapper('FL_PeakEM')
    PeakEMUnit = StringParamAttrWrapper('FL_PeakEMUnit')
    PeakFlux = StringParamAttrWrapper('FL_PeakFlux')
    PeakFluxUnit = StringParamAttrWrapper('FL_PeakFluxUnit')
    PeakTemp = StringParamAttrWrapper('FL_PeakTemp')
    PeakTempUnit = StringParamAttrWrapper('FL_PeakTempUnit')
    def __init__(self):
            ListAttr.__init__(self, "event_type", 'fl')

LP = ListAttr("event_type", 'lp')

OS = ListAttr("event_type", 'os')

@apply
class SS(ListAttr):
    SpinRate = StringParamAttrWrapper('SS_SpinRate')
    SpinRateUnit = StringParamAttrWrapper('SS_SpinRateUnit')
    def __init__(self):
            ListAttr.__init__(self, "event_type", 'ss')

@apply
class EF(ListAttr):
    AspectRatio = StringParamAttrWrapper('EF_AspectRatio')
    AxisLength = StringParamAttrWrapper('EF_AxisLength')
    AxisOrientation = StringParamAttrWrapper('EF_AxisOrientation')
    AxisOrientationUnit = StringParamAttrWrapper('EF_AxisOrientationUnit')
    FluxUnit = StringParamAttrWrapper('EF_FluxUnit')
    LengthUnit = StringParamAttrWrapper('EF_LengthUnit')
    NegEquivRadius = StringParamAttrWrapper('EF_NegEquivRadius')
    NegPeakFluxOnsetRate = StringParamAttrWrapper('EF_NegPeakFluxOnsetRate')
    OnsetRateUnit = StringParamAttrWrapper('EF_OnsetRateUnit')
    PosEquivRadius = StringParamAttrWrapper('EF_PosEquivRadius')
    PosPeakFluxOnsetRate = StringParamAttrWrapper('EF_PosPeakFluxOnsetRate')
    ProximityRatio = StringParamAttrWrapper('EF_ProximityRatio')
    SumNegSignedFlux = StringParamAttrWrapper('EF_SumNegSignedFlux')
    SumPosSignedFlux = StringParamAttrWrapper('EF_SumPosSignedFlux')
    def __init__(self):
            ListAttr.__init__(self, "event_type", 'ef')

CJ = ListAttr("event_type", 'cj')

PG = ListAttr("event_type", 'pg')

OT = ListAttr("event_type", 'ot')

NR = ListAttr("event_type", 'nr')

@apply
class SG(ListAttr):
    AspectRatio = StringParamAttrWrapper('SG_AspectRatio')
    Chirality = StringParamAttrWrapper('SG_Chirality')
    MeanContrast = StringParamAttrWrapper('SG_MeanContrast')
    Orientation = StringParamAttrWrapper('SG_Orientation')
    PeakContrast = StringParamAttrWrapper('SG_PeakContrast')
    Shape = StringParamAttrWrapper('SG_Shape')
    def __init__(self):
            ListAttr.__init__(self, "event_type", 'sg')

SP = ListAttr("event_type", 'sp')

CR = ListAttr("event_type", 'cr')

@apply
class CC(ListAttr):
    AxisUnit = StringParamAttrWrapper('CC_AxisUnit')
    MajorAxis = StringParamAttrWrapper('CC_MajorAxis')
    MinorAxis = StringParamAttrWrapper('CC_MinorAxis')
    TiltAngleMajorFromRadial = StringParamAttrWrapper('CC_TiltAngleMajorFromRadial')
    TiltAngleUnit = StringParamAttrWrapper('CC_TiltAngleUnit')
    def __init__(self):
            ListAttr.__init__(self, "event_type", 'cc')

ER = ListAttr("event_type", 'er')

@apply
class TO(ListAttr):
    Shape = StringParamAttrWrapper('TO_Shape')
    def __init__(self):
            ListAttr.__init__(self, "event_type", 'to')

class Misc(object):
    Area_AtDiskCenter = StringParamAttrWrapper('Area_AtDiskCenter')
    Area_AtDiskCenterUncert = StringParamAttrWrapper('Area_AtDiskCenterUncert')
    Area_Raw = StringParamAttrWrapper('Area_Raw')
    Area_Uncert = StringParamAttrWrapper('Area_Uncert')
    Area_Unit = StringParamAttrWrapper('Area_Unit')
    BoundBox_C1LL = StringParamAttrWrapper('BoundBox_C1LL')
    BoundBox_C1UR = StringParamAttrWrapper('BoundBox_C1UR')
    BoundBox_C2LL = StringParamAttrWrapper('BoundBox_C2LL')
    BoundBox_C2UR = StringParamAttrWrapper('BoundBox_C2UR')
    Bound_CCNsteps = StringParamAttrWrapper('Bound_CCNsteps')
    Bound_CCStartC1 = StringParamAttrWrapper('Bound_CCStartC1')
    Bound_CCStartC2 = StringParamAttrWrapper('Bound_CCStartC2')
    CME_Accel = StringParamAttrWrapper('CME_Accel')
    CME_AccelUncert = StringParamAttrWrapper('CME_AccelUncert')
    CME_AccelUnit = StringParamAttrWrapper('CME_AccelUnit')
    CME_AngularWidth = StringParamAttrWrapper('CME_AngularWidth')
    CME_AngularWidthUnit = StringParamAttrWrapper('CME_AngularWidthUnit')
    CME_Mass = StringParamAttrWrapper('CME_Mass')
    CME_MassUncert = StringParamAttrWrapper('CME_MassUncert')
    CME_MassUnit = StringParamAttrWrapper('CME_MassUnit')
    CME_RadialLinVel = StringParamAttrWrapper('CME_RadialLinVel')
    CME_RadialLinVelMax = StringParamAttrWrapper('CME_RadialLinVelMax')
    CME_RadialLinVelMin = StringParamAttrWrapper('CME_RadialLinVelMin')
    CME_RadialLinVelStddev = StringParamAttrWrapper('CME_RadialLinVelStddev')
    CME_RadialLinVelUncert = StringParamAttrWrapper('CME_RadialLinVelUncert')
    CME_RadialLinVelUnit = StringParamAttrWrapper('CME_RadialLinVelUnit')
    Event_C1Error = StringParamAttrWrapper('Event_C1Error')
    Event_C2Error = StringParamAttrWrapper('Event_C2Error')
    Event_ClippedSpatial = StringParamAttrWrapper('Event_ClippedSpatial')
    Event_ClippedTemporal = StringParamAttrWrapper('Event_ClippedTemporal')
    Event_Coord1 = StringParamAttrWrapper('Event_Coord1')
    Event_Coord2 = StringParamAttrWrapper('Event_Coord2')
    Event_Coord3 = StringParamAttrWrapper('Event_Coord3')
    Event_CoordSys = StringParamAttrWrapper('Event_CoordSys')
    Event_CoordUnit = StringParamAttrWrapper('Event_CoordUnit')
    Event_MapURL = StringParamAttrWrapper('Event_MapURL')
    Event_MaskURL = StringParamAttrWrapper('Event_MaskURL')
    Event_Npixels = StringParamAttrWrapper('Event_Npixels')
    Event_PixelUnit = StringParamAttrWrapper('Event_PixelUnit')
    Event_Probability = StringParamAttrWrapper('Event_Probability')
    Event_TestFlag = StringParamAttrWrapper('Event_TestFlag')
    Event_Type = StringParamAttrWrapper('Event_Type')
    FRM_Contact = StringParamAttrWrapper('FRM_Contact')
    FRM_HumanFlag = StringParamAttrWrapper('FRM_HumanFlag')
    FRM_Identifier = StringParamAttrWrapper('FRM_Identifier')
    FRM_Institute = StringParamAttrWrapper('FRM_Institute')
    FRM_Name = StringParamAttrWrapper('FRM_Name')
    FRM_ParamSet = StringParamAttrWrapper('FRM_ParamSet')
    FRM_SpecificID = StringParamAttrWrapper('FRM_SpecificID')
    FRM_URL = StringParamAttrWrapper('FRM_URL')
    FRM_VersionNumber = StringParamAttrWrapper('FRM_VersionNumber')
    FreqMaxRange = StringParamAttrWrapper('FreqMaxRange')
    FreqMinRange = StringParamAttrWrapper('FreqMinRange')
    FreqPeakPower = StringParamAttrWrapper('FreqPeakPower')
    FreqUnit = StringParamAttrWrapper('FreqUnit')
    IntensMaxAmpl = StringParamAttrWrapper('IntensMaxAmpl')
    IntensMinAmpl = StringParamAttrWrapper('IntensMinAmpl')
    IntensUnit = StringParamAttrWrapper('IntensUnit')
    KB_Archivist = StringParamAttrWrapper('KB_Archivist')
    MaxMagFieldStrength = StringParamAttrWrapper('MaxMagFieldStrength')
    MaxMagFieldStrengthUnit = StringParamAttrWrapper('MaxMagFieldStrengthUnit')
    OBS_ChannelID = StringParamAttrWrapper('OBS_ChannelID')
    OBS_DataPrepURL = StringParamAttrWrapper('OBS_DataPrepURL')
    OBS_FirstProcessingDate = StringParamAttrWrapper('OBS_FirstProcessingDate')
    OBS_IncludesNRT = StringParamAttrWrapper('OBS_IncludesNRT')
    OBS_Instrument = StringParamAttrWrapper('OBS_Instrument')
    OBS_LastProcessingDate = StringParamAttrWrapper('OBS_LastProcessingDate')
    OBS_LevelNum = StringParamAttrWrapper('OBS_LevelNum')
    OBS_MeanWavel = StringParamAttrWrapper('OBS_MeanWavel')
    OBS_Observatory = StringParamAttrWrapper('OBS_Observatory')
    OBS_Title = StringParamAttrWrapper('OBS_Title')
    OBS_WavelUnit = StringParamAttrWrapper('OBS_WavelUnit')
    OscillNPeriods = StringParamAttrWrapper('OscillNPeriods')
    OscillNPeriodsUncert = StringParamAttrWrapper('OscillNPeriodsUncert')
    Outflow_Length = StringParamAttrWrapper('Outflow_Length')
    Outflow_LengthUnit = StringParamAttrWrapper('Outflow_LengthUnit')
    Outflow_OpeningAngle = StringParamAttrWrapper('Outflow_OpeningAngle')
    Outflow_Speed = StringParamAttrWrapper('Outflow_Speed')
    Outflow_SpeedUnit = StringParamAttrWrapper('Outflow_SpeedUnit')
    Outflow_TransSpeed = StringParamAttrWrapper('Outflow_TransSpeed')
    Outflow_Width = StringParamAttrWrapper('Outflow_Width')
    Outflow_WidthUnit = StringParamAttrWrapper('Outflow_WidthUnit')
    PeakPower = StringParamAttrWrapper('PeakPower')
    PeakPowerUnit = StringParamAttrWrapper('PeakPowerUnit')
    RasterScanType = StringParamAttrWrapper('RasterScanType')
    Skel_Curvature = StringParamAttrWrapper('Skel_Curvature')
    Skel_Nsteps = StringParamAttrWrapper('Skel_Nsteps')
    Skel_StartC1 = StringParamAttrWrapper('Skel_StartC1')
    Skel_StartC2 = StringParamAttrWrapper('Skel_StartC2')
    VelocMaxAmpl = StringParamAttrWrapper('VelocMaxAmpl')
    VelocMaxPower = StringParamAttrWrapper('VelocMaxPower')
    VelocMaxPowerUncert = StringParamAttrWrapper('VelocMaxPowerUncert')
    VelocMinAmpl = StringParamAttrWrapper('VelocMinAmpl')
    VelocUnit = StringParamAttrWrapper('VelocUnit')
    WaveDisplMaxAmpl = StringParamAttrWrapper('WaveDisplMaxAmpl')
    WaveDisplMinAmpl = StringParamAttrWrapper('WaveDisplMinAmpl')
    WaveDisplUnit = StringParamAttrWrapper('WaveDisplUnit')
    WavelMaxPower = StringParamAttrWrapper('WavelMaxPower')
    WavelMaxPowerUncert = StringParamAttrWrapper('WavelMaxPowerUncert')
    WavelMaxRange = StringParamAttrWrapper('WavelMaxRange')
    WavelMinRange = StringParamAttrWrapper('WavelMinRange')
    WavelUnit = StringParamAttrWrapper('WavelUnit')
