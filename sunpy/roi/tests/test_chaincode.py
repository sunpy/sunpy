# Author: David PS <dps.helio -?- gmail.com>

import unittest
from sunpy.roi.chaincode import Chaincode
import numpy as np

class CCTests(unittest.TestCase):

    def testEnds(self):
        cc = Chaincode([0, 0], "2460") # Can I test more than one path? How?
        end = [0, 0]
        self.failUnless(cc.matchend(end))

    def testEndsFalse(self):
        cc = Chaincode([0, 0], "24460")
        end = [0, 0]
        self.failIf(cc.matchend(end))

    def testSecondCoordinate(self):
        cc = Chaincode([0, 0], "0023")
        second = [-2, 0]
        self.failUnless(cc.matchany(second, 2))

    def testSecondCoordinateFails(self):
        cc = Chaincode([1, 0], "0023")
        second = [-2, 0]
        self.failIf(cc.matchany(second, 2))

    def testScaleSecond(self):
        cc = Chaincode([0, 0], "0723", xdelta=0.5, ydelta=0.5)
        second = [-1, 0.5]
        self.failUnless(cc.matchany(second, 2))

    def testScaleEnd(self):
        cc = Chaincode([1.2, 3],"0723", xdelta=2.629, ydelta=2.629)
        end = [-1.429, 0.371]
        self.failUnless(cc.matchany(end, -1))

    def testnparray(self):
        # Let's test that the shape of the array matches the expected
        # To do so we need to use np.array, instead of lists.
        cc = Chaincode([0, 0], "2460")
        shape = (2,5)
        self.failUnless(cc.coordinates.shape == shape)

    def testBoundingBox(self): #needs of np.array... I think
        cc = Chaincode([0, 0], "00033344")
        boundingbox = [[-3, 2], [-3, 0]] # [[x0,x1],[y0,y1]] (like cc)
        self.failUnless(np.all(cc.BoundingBox() == np.array(boundingbox)))

    def testBoundingBoxFalse(self):
        cc = Chaincode([0, 0], "002")
        boundingbox = [[-1, 0], [-1, 0]]
        self.failIf(np.all(cc.BoundingBox() != np.array(boundingbox)))

    def testSubBoundingBoxX(self):
        cc = Chaincode([0, 0], "44464660012075602223")
        self.failUnless(cc.subBoundingBox(xedge=[0.1, 2]) == [0, 3])

    def testSubBoundingBoxY(self):
        cc = Chaincode([0, 0], "44464660012075602223")
        self.failUnless(cc.subBoundingBox(yedge=[-1, 0.5]) == [0, 3])



def main():
    unittest.main()

if __name__ == '__main__':
    main()
