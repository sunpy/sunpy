# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

# pylint: disable=C0103,R0903

from urllib2 import urlopen
from urllib import urlencode
from datetime import datetime
from functools import partial

from sunpy.net import attr

DEFAULT_URL = 'http://www.lmsal.com/hek/her'


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
    dct['event_starttime'] = root.start.strftime('%Y-%m-%dT%H:%M:%S')
    dct['event_endtime'] = root.end.strftime('%Y-%m-%dT%H:%M:%S')
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

EVENTS = [
    'AR', 'CE', 'CD', 'CH', 'CW', 'FI', 'FE', 'FA', 'FL', 'LP', 'OS', 'SS',
    'EF', 'CJ', 'PG', 'OT', 'NR', 'SG', 'SP', 'CR', 'CC', 'ER', 'TO'
]

    
class HEKClient(object):
    # FIXME: Types!
    fields = {
        'AR_CompactnessCls': StringParamAttrWrapper,
        'AR_IntensKurt': StringParamAttrWrapper,
        'AR_IntensMax': StringParamAttrWrapper,
        'AR_IntensMean': StringParamAttrWrapper,
        'AR_IntensMin': StringParamAttrWrapper,
        'AR_IntensSkew': StringParamAttrWrapper,
        'AR_IntensTotal': StringParamAttrWrapper,
        'AR_IntensUnit': StringParamAttrWrapper,
        'AR_IntensVar': StringParamAttrWrapper,
        'AR_McIntoshCls': StringParamAttrWrapper,
        'AR_MtWilsonCls': StringParamAttrWrapper,
        'AR_NOAANum': StringParamAttrWrapper,
        'AR_NOAAclass': StringParamAttrWrapper,
        'AR_NumSpots': StringParamAttrWrapper,
        'AR_PenumbraCls': StringParamAttrWrapper,
        'AR_Polarity': StringParamAttrWrapper,
        'AR_SpotAreaRaw': StringParamAttrWrapper,
        'AR_SpotAreaRawUncert': StringParamAttrWrapper,
        'AR_SpotAreaRawUnit': StringParamAttrWrapper,
        'AR_SpotAreaRepr': StringParamAttrWrapper,
        'AR_SpotAreaReprUncert': StringParamAttrWrapper,
        'AR_SpotAreaReprUnit': StringParamAttrWrapper,
        'AR_ZurichCls': StringParamAttrWrapper,
        'Area_AtDiskCenter': StringParamAttrWrapper,
        'Area_AtDiskCenterUncert': StringParamAttrWrapper,
        'Area_Raw': StringParamAttrWrapper,
        'Area_Uncert': StringParamAttrWrapper,
        'Area_Unit': StringParamAttrWrapper,
        'BoundBox_C1LL': StringParamAttrWrapper,
        'BoundBox_C1UR': StringParamAttrWrapper,
        'BoundBox_C2LL': StringParamAttrWrapper,
        'BoundBox_C2UR': StringParamAttrWrapper,
        'Bound_CCNsteps': StringParamAttrWrapper,
        'Bound_CCStartC1': StringParamAttrWrapper,
        'Bound_CCStartC2': StringParamAttrWrapper,
        'CC_AxisUnit': StringParamAttrWrapper,
        'CC_MajorAxis': StringParamAttrWrapper,
        'CC_MinorAxis': StringParamAttrWrapper,
        'CC_TiltAngleMajorFromRadial': StringParamAttrWrapper,
        'CC_TiltAngleUnit': StringParamAttrWrapper,
        'CD_Area': StringParamAttrWrapper,
        'CD_AreaUncert': StringParamAttrWrapper,
        'CD_AreaUnit': StringParamAttrWrapper,
        'CD_Mass': StringParamAttrWrapper,
        'CD_MassUncert': StringParamAttrWrapper,
        'CD_MassUnit': StringParamAttrWrapper,
        'CD_Volume': StringParamAttrWrapper,
        'CD_VolumeUncert': StringParamAttrWrapper,
        'CD_VolumeUnit': StringParamAttrWrapper,
        'CME_Accel': StringParamAttrWrapper,
        'CME_AccelUncert': StringParamAttrWrapper,
        'CME_AccelUnit': StringParamAttrWrapper,
        'CME_AngularWidth': StringParamAttrWrapper,
        'CME_AngularWidthUnit': StringParamAttrWrapper,
        'CME_Mass': StringParamAttrWrapper,
        'CME_MassUncert': StringParamAttrWrapper,
        'CME_MassUnit': StringParamAttrWrapper,
        'CME_RadialLinVel': StringParamAttrWrapper,
        'CME_RadialLinVelMax': StringParamAttrWrapper,
        'CME_RadialLinVelMin': StringParamAttrWrapper,
        'CME_RadialLinVelStddev': StringParamAttrWrapper,
        'CME_RadialLinVelUncert': StringParamAttrWrapper,
        'CME_RadialLinVelUnit': StringParamAttrWrapper,
        'EF_AspectRatio': StringParamAttrWrapper,
        'EF_AxisLength': StringParamAttrWrapper,
        'EF_AxisOrientation': StringParamAttrWrapper,
        'EF_AxisOrientationUnit': StringParamAttrWrapper,
        'EF_FluxUnit': StringParamAttrWrapper,
        'EF_LengthUnit': StringParamAttrWrapper,
        'EF_NegEquivRadius': StringParamAttrWrapper,
        'EF_NegPeakFluxOnsetRate': StringParamAttrWrapper,
        'EF_OnsetRateUnit': StringParamAttrWrapper,
        'EF_PosEquivRadius': StringParamAttrWrapper,
        'EF_PosPeakFluxOnsetRate': StringParamAttrWrapper,
        'EF_ProximityRatio': StringParamAttrWrapper,
        'EF_SumNegSignedFlux': StringParamAttrWrapper,
        'EF_SumPosSignedFlux': StringParamAttrWrapper,
        'Event_C1Error': StringParamAttrWrapper,
        'Event_C2Error': StringParamAttrWrapper,
        'Event_ClippedSpatial': StringParamAttrWrapper,
        'Event_ClippedTemporal': StringParamAttrWrapper,
        'Event_Coord1': StringParamAttrWrapper,
        'Event_Coord2': StringParamAttrWrapper,
        'Event_Coord3': StringParamAttrWrapper,
        'Event_CoordSys': StringParamAttrWrapper,
        'Event_CoordUnit': StringParamAttrWrapper,
        'Event_MapURL': StringParamAttrWrapper,
        'Event_MaskURL': StringParamAttrWrapper,
        'Event_Npixels': StringParamAttrWrapper,
        'Event_PixelUnit': StringParamAttrWrapper,
        'Event_Probability': StringParamAttrWrapper,
        'Event_TestFlag': StringParamAttrWrapper,
        'Event_Type': StringParamAttrWrapper,
        'FI_BarbsL': StringParamAttrWrapper,
        'FI_BarbsR': StringParamAttrWrapper,
        'FI_BarbsTot': StringParamAttrWrapper,
        'FI_Chirality': StringParamAttrWrapper,
        'FI_Length': StringParamAttrWrapper,
        'FI_LengthUnit': StringParamAttrWrapper,
        'FI_Tilt': StringParamAttrWrapper,
        'FL_EFoldTime': StringParamAttrWrapper,
        'FL_EFoldTimeUnit': StringParamAttrWrapper,
        'FL_Fluence': StringParamAttrWrapper,
        'FL_FluenceUnit': StringParamAttrWrapper,
        'FL_GOESCls': StringParamAttrWrapper,
        'FL_PeakEM': StringParamAttrWrapper,
        'FL_PeakEMUnit': StringParamAttrWrapper,
        'FL_PeakFlux': StringParamAttrWrapper,
        'FL_PeakFluxUnit': StringParamAttrWrapper,
        'FL_PeakTemp': StringParamAttrWrapper,
        'FL_PeakTempUnit': StringParamAttrWrapper,
        'FRM_Contact': StringParamAttrWrapper,
        'FRM_HumanFlag': StringParamAttrWrapper,
        'FRM_Identifier': StringParamAttrWrapper,
        'FRM_Institute': StringParamAttrWrapper,
        'FRM_Name': StringParamAttrWrapper,
        'FRM_ParamSet': StringParamAttrWrapper,
        'FRM_SpecificID': StringParamAttrWrapper,
        'FRM_URL': StringParamAttrWrapper,
        'FRM_VersionNumber': StringParamAttrWrapper,
        'FreqMaxRange': StringParamAttrWrapper,
        'FreqMinRange': StringParamAttrWrapper,
        'FreqPeakPower': StringParamAttrWrapper,
        'FreqUnit': StringParamAttrWrapper,
        'IntensMaxAmpl': StringParamAttrWrapper,
        'IntensMinAmpl': StringParamAttrWrapper,
        'IntensUnit': StringParamAttrWrapper,
        'KB_Archivist': StringParamAttrWrapper,
        'MaxMagFieldStrength': StringParamAttrWrapper,
        'MaxMagFieldStrengthUnit': StringParamAttrWrapper,
        'OBS_ChannelID': StringParamAttrWrapper,
        'OBS_DataPrepURL': StringParamAttrWrapper,
        'OBS_FirstProcessingDate': StringParamAttrWrapper,
        'OBS_IncludesNRT': StringParamAttrWrapper,
        'OBS_Instrument': StringParamAttrWrapper,
        'OBS_LastProcessingDate': StringParamAttrWrapper,
        'OBS_LevelNum': StringParamAttrWrapper,
        'OBS_MeanWavel': StringParamAttrWrapper,
        'OBS_Observatory': StringParamAttrWrapper,
        'OBS_Title': StringParamAttrWrapper,
        'OBS_WavelUnit': StringParamAttrWrapper,
        'OscillNPeriods': StringParamAttrWrapper,
        'OscillNPeriodsUncert': StringParamAttrWrapper,
        'Outflow_Length': StringParamAttrWrapper,
        'Outflow_LengthUnit': StringParamAttrWrapper,
        'Outflow_OpeningAngle': StringParamAttrWrapper,
        'Outflow_Speed': StringParamAttrWrapper,
        'Outflow_SpeedUnit': StringParamAttrWrapper,
        'Outflow_TransSpeed': StringParamAttrWrapper,
        'Outflow_Width': StringParamAttrWrapper,
        'Outflow_WidthUnit': StringParamAttrWrapper,
        'PeakPower': StringParamAttrWrapper,
        'PeakPowerUnit': StringParamAttrWrapper,
        'RasterScanType': StringParamAttrWrapper,
        'SG_AspectRatio': StringParamAttrWrapper,
        'SG_Chirality': StringParamAttrWrapper,
        'SG_MeanContrast': StringParamAttrWrapper,
        'SG_Orientation': StringParamAttrWrapper,
        'SG_PeakContrast': StringParamAttrWrapper,
        'SG_Shape': StringParamAttrWrapper,
        'SS_SpinRate': StringParamAttrWrapper,
        'SS_SpinRateUnit': StringParamAttrWrapper,
        'Skel_Curvature': StringParamAttrWrapper,
        'Skel_Nsteps': StringParamAttrWrapper,
        'Skel_StartC1': StringParamAttrWrapper,
        'Skel_StartC2': StringParamAttrWrapper,
        'TO_Shape': StringParamAttrWrapper,
        'VelocMaxAmpl': StringParamAttrWrapper,
        'VelocMaxPower': StringParamAttrWrapper,
        'VelocMaxPowerUncert': StringParamAttrWrapper,
        'VelocMinAmpl': StringParamAttrWrapper,
        'VelocUnit': StringParamAttrWrapper,
        'WaveDisplMaxAmpl': StringParamAttrWrapper,
        'WaveDisplMinAmpl': StringParamAttrWrapper,
        'WaveDisplUnit': StringParamAttrWrapper,
        'WavelMaxPower': StringParamAttrWrapper,
        'WavelMaxPowerUncert': StringParamAttrWrapper,
        'WavelMaxRange': StringParamAttrWrapper,
        'WavelMinRange': StringParamAttrWrapper,
        'WavelUnit': StringParamAttrWrapper
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
    walker.apply(SpartialRegion(), {}, default)
    
    def __init__(self, url=DEFAULT_URL):
        self.url = url
    
    def _download(self, data):
        page = 1        
        results = []
        
        while True:
            data['page'] = page
            result = json.load(urlopen(self.url, urlencode(data)))
            results.extend(result['result'])
            
            if not result['overmax']:
                return results
            page += 1
    
    def query(self, *query):
        query = attr.and_(*query)
        
        data = walker.create(query, {})
        ndata = []
        for elem in data:
            new = self.default.copy()
            new.update(elem)
            ndata.append(new)
        
        if len(ndata) == 1:
            return self._download(ndata[0])
        else:
            return self.merge(self._download(data) for data in ndata)
    
    def merge(self, responses):
        return responses
    
    def __getitem__(self, name):
        return self.fields[name](name)


if __name__ == '__main__':
    import json
    import pprint
    c = HEKClient()
    print len(c.query(
        Time.dt((2010, 1, 1), (2010, 1, 23, 23)),
        c['AR'], c['FL']
    ))
