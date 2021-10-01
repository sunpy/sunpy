import numpy as np

from sunpy.net.helio.chaincode import Chaincode


def test_ends():
    cc = Chaincode([0, 0], "2460")
    end = [0, 0]
    assert cc.matchend(end)


def test_endsfalse():
    cc = Chaincode([0, 0], "24460")
    end = [0, 0]
    assert not cc.matchend(end)


def test_secondcoordinate():
    cc = Chaincode([0, 0], "0023")
    second = [-2, 0]
    assert cc.matchany(second, 2)


def test_secondcoordinatefails():
    cc = Chaincode([1, 0], "0023")
    second = [-2, 0]
    assert not cc.matchany(second, 2)


def test_scalesecond():
    cc = Chaincode([0, 0], "0723", xdelta=0.5, ydelta=0.5)
    second = [-1, 0.5]
    assert cc.matchany(second, 2)


def test_scaleend():
    cc = Chaincode([1.2, 3], "0723", xdelta=2.629, ydelta=2.629)
    end = [-1.429, 0.371]
    assert cc.matchany(end, -1)


def test_array():
    cc = Chaincode(np.array([0, 0]), "2460")
    shape = (2, 5)
    assert cc.coordinates.shape == shape


def test_boundingbox():
    cc = Chaincode([0, 0], "00033344")
    boundingbox = [[-3, 2], [-3, 0]]  # [[x0,x1],[y0,y1]] (like cc)
    assert np.all(cc.boundingbox() == np.array(boundingbox))


def test_boundingboxfalse():
    cc = Chaincode([0, 0], "002")
    boundingbox = [[-1, 0], [-1, 0]]
    assert not np.all(cc.boundingbox() != np.array(boundingbox))


def test_subboundingboxx():
    cc = Chaincode([0, 0], "44464660012075602223")
    assert cc.sub_boundingbox(xedge=[0.1, 2]) == [0, 3]


def test_subboundingboxy():
    cc = Chaincode([0, 0], "44464660012075602223")
    assert cc.sub_boundingbox(yedge=[-1, 0.5]) == [0, 3]
