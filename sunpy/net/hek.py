# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

import sys

from urllib2 import urlopen
from urllib import urlencode
from datetime import datetime
from functools import partial

from sunpy.net import attr

DEFAULT_URL = 'http://www.lmsal.com/hek/her'


class ParamAttr(attr.KeysAttr):
    def __init__(self, name, op, value):
        attr.KeysAttr.__init__(self, [(name, op)])
        self.name = name
        self.op = op
        self.value = value


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


class SpartialRegion(attr.Attr):
    def __init__(
        self, x1=-1200, y1=-1200, x2=1200, y2=1200, sys='helioprojective'):
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
        self.types = types
    
    def collides(self, other):
        return False
    
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return vars(self) == vars(other)
    
    def __hash__(self):
        return hash(tuple(vars(self).itervalues()))


walker = attr.AttrWalker()


@walker.add_applier(Contains)
def _a(walker, root, state, dct):
    dct['type'] = 'contains'
    if not Contains in state:
        state[Contains] = 1
    
    nid = state[Contains]
    for n, type_ in enumerate(root.types):
        dct['event_type%d' % (nid + n)] = type_
    state[Contains] += n
    return dct

@walker.add_creator(
    Time, SpartialRegion, ListAttr, ParamAttr, attr.AttrAnd, Contains)
def _c(walker, root, state):
    value = {}
    walker.apply(root, state, value)
    return [value]

@walker.add_applier(Time)
def _a(walker, root, state, dct):
    dct['event_starttime'] = root.start.strftime('%Y-%m-%dT%H:%M:%S')
    dct['event_endtime'] = root.end.strftime('%Y-%m-%dT%H:%M:%S')
    return dct

@walker.add_applier(SpartialRegion)
def _a(walker, root, state, dct):
    dct['x1'] = root.x1
    dct['y1'] = root.y1
    dct['x2'] = root.x2
    dct['y2'] = root.y2
    dct['event_coordsys'] = root.sys
    return dct

@walker.add_applier(ListAttr)
def _a(walker, root, state, dct):
    if root.key in dct:
        dct[root.key] += ',%s' % root.item
    else:
        dct[root.key] = root.item
    return dct

@walker.add_applier(EventType)
def _a(walker, root, state, dct):
    if dct.get('type', None) == 'contains':
        raise ValueError
    
    return walker.super_apply(super(EventType, root), state, dct)

@walker.add_applier(ParamAttr)
def _a(walker, root, state, dct):
    if not ParamAttr in state:
        state[ParamAttr] = 0
    
    nid = state[ParamAttr]
    dct['param%d' % nid] = root.name
    dct['op%d' % nid] = root.op
    dct['value%d' % nid] = root.value
    state[ParamAttr] += 1
    return dct

@walker.add_applier(attr.AttrAnd)
def _a(walker, root, state, dct):
    for attr in root.attrs:
        walker.apply(attr, state, dct)

@walker.add_creator(attr.AttrOr)
def _c(walker, root, state):
    blocks = []
    for attr in root.attrs:
        blocks.extend(walker.create(attr, state))
    return blocks

@walker.add_creator(attr.DummyAttr)
def _c(walker, root, state):
    return {}

@walker.add_applier(attr.DummyAttr)
def _a(walker, root, state, dct):
    pass


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
    def __init__(self, name):
        self.name = name
    
    def like(self, other):
        return ParamAttr(self.name, 'like', other)


class NumberParamAttrWrapper(ComparisonParamAttrWrapper):
    pass

EVENTS = [
    'AR', 'CE', 'CD', 'CH', 'CW', 'FI', 'FE', 'FA', 'FL', 'LP', 'OS', 'SS',
    'EF', 'CJ', 'PG', 'OT', 'NR', 'SG', 'SP', 'CR', 'CC', 'ER', 'TO'
]

    
class HEKClient(object):
    # FIXME: Types!
    fields = {
        'AR_CompactnessCls': StringParamWrapper,
        'AR_IntensKurt': StringParamWrapper,
        'AR_IntensMax': StringParamWrapper,
        'AR_IntensMean': StringParamWrapper,
        'AR_IntensMin': StringParamWrapper,
        'AR_IntensSkew': StringParamWrapper,
        'AR_IntensTotal': StringParamWrapper,
        'AR_IntensUnit': StringParamWrapper,
        'AR_IntensVar': StringParamWrapper,
        'AR_McIntoshCls': StringParamWrapper,
        'AR_MtWilsonCls': StringParamWrapper,
        'AR_NOAANum': StringParamWrapper,
        'AR_NOAAclass': StringParamWrapper,
        'AR_NumSpots': StringParamWrapper,
        'AR_PenumbraCls': StringParamWrapper,
        'AR_Polarity': StringParamWrapper,
        'AR_SpotAreaRaw': StringParamWrapper,
        'AR_SpotAreaRawUncert': StringParamWrapper,
        'AR_SpotAreaRawUnit': StringParamWrapper,
        'AR_SpotAreaRepr': StringParamWrapper,
        'AR_SpotAreaReprUncert': StringParamWrapper,
        'AR_SpotAreaReprUnit': StringParamWrapper,
        'AR_ZurichCls': StringParamWrapper,
        'Area_AtDiskCenter': StringParamWrapper,
        'Area_AtDiskCenterUncert': StringParamWrapper,
        'Area_Raw': StringParamWrapper,
        'Area_Uncert': StringParamWrapper,
        'Area_Unit': StringParamWrapper,
        'BoundBox_C1LL': StringParamWrapper,
        'BoundBox_C1UR': StringParamWrapper,
        'BoundBox_C2LL': StringParamWrapper,
        'BoundBox_C2UR': StringParamWrapper,
        'Bound_CCNsteps': StringParamWrapper,
        'Bound_CCStartC1': StringParamWrapper,
        'Bound_CCStartC2': StringParamWrapper,
        'CC_AxisUnit': StringParamWrapper,
        'CC_MajorAxis': StringParamWrapper,
        'CC_MinorAxis': StringParamWrapper,
        'CC_TiltAngleMajorFromRadial': StringParamWrapper,
        'CC_TiltAngleUnit': StringParamWrapper,
        'CD_Area': StringParamWrapper,
        'CD_AreaUncert': StringParamWrapper,
        'CD_AreaUnit': StringParamWrapper,
        'CD_Mass': StringParamWrapper,
        'CD_MassUncert': StringParamWrapper,
        'CD_MassUnit': StringParamWrapper,
        'CD_Volume': StringParamWrapper,
        'CD_VolumeUncert': StringParamWrapper,
        'CD_VolumeUnit': StringParamWrapper,
        'CME_Accel': StringParamWrapper,
        'CME_AccelUncert': StringParamWrapper,
        'CME_AccelUnit': StringParamWrapper,
        'CME_AngularWidth': StringParamWrapper,
        'CME_AngularWidthUnit': StringParamWrapper,
        'CME_Mass': StringParamWrapper,
        'CME_MassUncert': StringParamWrapper,
        'CME_MassUnit': StringParamWrapper,
        'CME_RadialLinVel': StringParamWrapper,
        'CME_RadialLinVelMax': StringParamWrapper,
        'CME_RadialLinVelMin': StringParamWrapper,
        'CME_RadialLinVelStddev': StringParamWrapper,
        'CME_RadialLinVelUncert': StringParamWrapper,
        'CME_RadialLinVelUnit': StringParamWrapper,
        'EF_AspectRatio': StringParamWrapper,
        'EF_AxisLength': StringParamWrapper,
        'EF_AxisOrientation': StringParamWrapper,
        'EF_AxisOrientationUnit': StringParamWrapper,
        'EF_FluxUnit': StringParamWrapper,
        'EF_LengthUnit': StringParamWrapper,
        'EF_NegEquivRadius': StringParamWrapper,
        'EF_NegPeakFluxOnsetRate': StringParamWrapper,
        'EF_OnsetRateUnit': StringParamWrapper,
        'EF_PosEquivRadius': StringParamWrapper,
        'EF_PosPeakFluxOnsetRate': StringParamWrapper,
        'EF_ProximityRatio': StringParamWrapper,
        'EF_SumNegSignedFlux': StringParamWrapper,
        'EF_SumPosSignedFlux': StringParamWrapper,
        'Event_C1Error': StringParamWrapper,
        'Event_C2Error': StringParamWrapper,
        'Event_ClippedSpatial': StringParamWrapper,
        'Event_ClippedTemporal': StringParamWrapper,
        'Event_Coord1': StringParamWrapper,
        'Event_Coord2': StringParamWrapper,
        'Event_Coord3': StringParamWrapper,
        'Event_CoordSys': StringParamWrapper,
        'Event_CoordUnit': StringParamWrapper,
        'Event_MapURL': StringParamWrapper,
        'Event_MaskURL': StringParamWrapper,
        'Event_Npixels': StringParamWrapper,
        'Event_PixelUnit': StringParamWrapper,
        'Event_Probability': StringParamWrapper,
        'Event_TestFlag': StringParamWrapper,
        'Event_Type': StringParamWrapper,
        'FI_BarbsL': StringParamWrapper,
        'FI_BarbsR': StringParamWrapper,
        'FI_BarbsTot': StringParamWrapper,
        'FI_Chirality': StringParamWrapper,
        'FI_Length': StringParamWrapper,
        'FI_LengthUnit': StringParamWrapper,
        'FI_Tilt': StringParamWrapper,
        'FL_EFoldTime': StringParamWrapper,
        'FL_EFoldTimeUnit': StringParamWrapper,
        'FL_Fluence': StringParamWrapper,
        'FL_FluenceUnit': StringParamWrapper,
        'FL_GOESCls': StringParamWrapper,
        'FL_PeakEM': StringParamWrapper,
        'FL_PeakEMUnit': StringParamWrapper,
        'FL_PeakFlux': StringParamWrapper,
        'FL_PeakFluxUnit': StringParamWrapper,
        'FL_PeakTemp': StringParamWrapper,
        'FL_PeakTempUnit': StringParamWrapper,
        'FRM_Contact': StringParamWrapper,
        'FRM_HumanFlag': StringParamWrapper,
        'FRM_Identifier': StringParamWrapper,
        'FRM_Institute': StringParamWrapper,
        'FRM_Name': StringParamWrapper,
        'FRM_ParamSet': StringParamWrapper,
        'FRM_SpecificID': StringParamWrapper,
        'FRM_URL': StringParamWrapper,
        'FRM_VersionNumber': StringParamWrapper,
        'FreqMaxRange': StringParamWrapper,
        'FreqMinRange': StringParamWrapper,
        'FreqPeakPower': StringParamWrapper,
        'FreqUnit': StringParamWrapper,
        'IntensMaxAmpl': StringParamWrapper,
        'IntensMinAmpl': StringParamWrapper,
        'IntensUnit': StringParamWrapper,
        'KB_Archivist': StringParamWrapper,
        'MaxMagFieldStrength': StringParamWrapper,
        'MaxMagFieldStrengthUnit': StringParamWrapper,
        'OBS_ChannelID': StringParamWrapper,
        'OBS_DataPrepURL': StringParamWrapper,
        'OBS_FirstProcessingDate': StringParamWrapper,
        'OBS_IncludesNRT': StringParamWrapper,
        'OBS_Instrument': StringParamWrapper,
        'OBS_LastProcessingDate': StringParamWrapper,
        'OBS_LevelNum': StringParamWrapper,
        'OBS_MeanWavel': StringParamWrapper,
        'OBS_Observatory': StringParamWrapper,
        'OBS_Title': StringParamWrapper,
        'OBS_WavelUnit': StringParamWrapper,
        'OscillNPeriods': StringParamWrapper,
        'OscillNPeriodsUncert': StringParamWrapper,
        'Outflow_Length': StringParamWrapper,
        'Outflow_LengthUnit': StringParamWrapper,
        'Outflow_OpeningAngle': StringParamWrapper,
        'Outflow_Speed': StringParamWrapper,
        'Outflow_SpeedUnit': StringParamWrapper,
        'Outflow_TransSpeed': StringParamWrapper,
        'Outflow_Width': StringParamWrapper,
        'Outflow_WidthUnit': StringParamWrapper,
        'PeakPower': StringParamWrapper,
        'PeakPowerUnit': StringParamWrapper,
        'RasterScanType': StringParamWrapper,
        'SG_AspectRatio': StringParamWrapper,
        'SG_Chirality': StringParamWrapper,
        'SG_MeanContrast': StringParamWrapper,
        'SG_Orientation': StringParamWrapper,
        'SG_PeakContrast': StringParamWrapper,
        'SG_Shape': StringParamWrapper,
        'SS_SpinRate': StringParamWrapper,
        'SS_SpinRateUnit': StringParamWrapper,
        'Skel_Curvature': StringParamWrapper,
        'Skel_Nsteps': StringParamWrapper,
        'Skel_StartC1': StringParamWrapper,
        'Skel_StartC2': StringParamWrapper,
        'TO_Shape': StringParamWrapper,
        'VelocMaxAmpl': StringParamWrapper,
        'VelocMaxPower': StringParamWrapper,
        'VelocMaxPowerUncert': StringParamWrapper,
        'VelocMinAmpl': StringParamWrapper,
        'VelocUnit': StringParamWrapper,
        'WaveDisplMaxAmpl': StringParamWrapper,
        'WaveDisplMinAmpl': StringParamWrapper,
        'WaveDisplUnit': StringParamWrapper,
        'WavelMaxPower': StringParamWrapper,
        'WavelMaxPowerUncert': StringParamWrapper,
        'WavelMaxRange': StringParamWrapper,
        'WavelMinRange': StringParamWrapper,
        'WavelUnit': StringParamWrapper
    }
    
    for elem in EVENTS:
        fields[elem] = partial(ListAttr, 'event_type')
    
    default = {
        'cosec': '2',
        'cmd': 'search',
        'type': 'column',
        'event_type': '**',
    }
    # Default to full disk.
    walker.apply(SpartialRegion(), [], default)
    
    def __init__(self, url=DEFAULT_URL):
        self.url = url
    
    def query(self, *query):
        if len(query) > 1:
            query = attr.and_(*query)
        
        data = walker.create(query, {})
        ndata = []
        for elem in data:
            new = self.default.copy()
            new.update(elem)
            ndata.append(new)
        
        return urlopen(self.url, urlencode(ndata[0]))
    
    def __getitem__(self, name):
        return self.fields[name](name)


if __name__ == '__main__':
    import json
    import pprint
    c = HEKClient()
    pprint.pprint(json.load(c.query(
        Time.dt((2010, 1, 1), (2010, 1, 1, 1)),
        c['AR'],
    )))
