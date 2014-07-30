
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
