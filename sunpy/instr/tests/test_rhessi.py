
import datetime
import pytest
from sunpy.instr import rhessi
from sunpy.time import parse_time


@pytest.mark.online
def test_rhessi_flare_list():
    #create the flare list object
    flares = rhessi.RHESSIFlareList()
    #test some entries of the flare list
    assert type(flares.flare_list) == list
    assert flares.flare_list[0]['Flare ID'] == '2021213'
    assert parse_time(flares.flare_list[0]['Start date']) == datetime.datetime(2002, 2, 12, 0, 0)
    assert flares.flare_list[0]['Dur (s)'] == 712.0
    assert flares.flare_list[0]['Peak c/s'] == 136.0
    assert flares.flare_list[0]['Total Counts'] == 167304.0
    assert flares.flare_list[0]['Energy (keV)'] == '12-25'
    assert flares.flare_list[0]['X pos (arcsec)'] == 592.0
    assert flares.flare_list[0]['Y pos (arcsec)'] == -358.0
    assert flares.flare_list[0]['Radial (arcsec)'] == 692.0
    assert flares.flare_list[0]['AR'] == '0'
    assert flares.flare_list[0]['Flags'] == ['A1', 'P1']


@pytest.mark.online
def test_rhessi_flare_list_by_date_range():
    flares = rhessi.RHESSIFlarelist()
    #extract a selection of flares from the RHESSI flare list using a date range
    subflares = flares.find_events_by_date_range('2011-02-01','2011-02-03')
    assert len(subflares) == 14
    assert type(subflares) == list
    assert subflares[0]['Flare ID'] == '11020101'
    assert parse_time(subflares[0]['Start date']) == datetime.datetime(2011, 2, 1, 0, 0)
    assert subflares[0]['Dur (s)'] == 588.0
    assert subflares[0]['Peak c/s'] == 64.0
    assert subflares[0]['Total Counts'] == 108672.0
    assert subflares[0]['Energy (keV)'] == '6-12'
    assert subflares[0]['X pos (arcsec)'] == -337.0
    assert subflares[0]['Y pos (arcsec)'] == -238.0
    assert subflares[0]['Radial (arcsec)'] == 413.0
    assert subflares[0]['AR'] == '1150'
    assert subflares[0]['Flags'] == ['A0', 'DF', 'P1', 'Q1']
