"""
This module implements a function to understand chaincodes from HELIO.
"""
import numpy as np

__all__ = ["Chaincode"]


class Chaincode(np.ndarray):
    """
    A tool to infer some information from chaincodes produced by HELIO Feature
    Catalogue or Heliophysics Events Knowledgebase.

    Parameters
    ----------
    origin : `numpy.ndarray`, `list`
        The 2 points of the origin of the chaincode.
    chaincode : `str`
        A list of the numbers (0-7) that indicate the path of the chaincode.
        0 moves horizontally to the left and the rest follows anticlockwise.
    xdelta : `float`, optional
        The scale to convert between pixels and flat coordinates. Defaults to 1.0.
    ydelta : `float`, optional
        The scale to convert between pixels and flat coordinates. Defaults to 1.0

    Returns
    -------
    `numpy.ndarray`
        An array containing all the x and y coordinates of the chaincode
        e.g., ``[[x0, x1, x2, ..., xn], [y0 ,y1, y2, ..., yn]]``.

    Examples
    --------
    >>> from sunpy.net.helio.chaincode import Chaincode
    >>> cc = Chaincode([-88, 812], "44464655567670006011212222324",
    ...     xdelta=2.629, ydelta=2.629)
    >>> fig = plt.figure()   # doctest: +SKIP
    >>> ax = fig.add_subplot(111)   # doctest: +SKIP
    >>> x,y = zip(cc.coordinates)   # doctest: +SKIP
    >>> ax.plot(x[0], y[0], 'go-')   # doctest: +SKIP
    >>> fig.show()   # doctest: +SKIP
    """
    def __new__(cls, origin, chaincode, **kwargs):
        if isinstance(origin, list):
            obj = np.asarray(origin).view(cls)
        elif isinstance(origin, np.ndarray):
            obj = origin.view(cls)
        else:
            raise TypeError('Invalid input')
        return obj

    def __init__(self, origin, chaincode, xdelta=1, ydelta=1):
        x_steps = [-1, -1, 0, 1, 1, 1, 0, -1]
        y_steps = [0, -1, -1, -1, 0, 1, 1, 1]
        self.coordinates = np.ndarray((2, len(chaincode) + 1))
        self.coordinates[:, 0] = origin
        if chaincode.isdigit():
            for index, step in enumerate(chaincode):
                self.coordinates[:, index + 1] = self.coordinates[:, index] + \
                    [[x_steps[int(step)] * xdelta,
                      y_steps[int(step)] * ydelta]]

    def matchend(self, end):
        return np.all(np.equal(self.coordinates[:, -1], np.asarray(end)))

    def matchany(self, coordinates, index):
        return np.all(
            np.allclose(self.coordinates[:, index], np.asarray(coordinates))
        )

    def boundingbox(self):
        """
        Extract the coordinates of the chaincode.

        Returns in the form of ``[[x0, x1], [y0, y1]]``.
        """
        bb = np.zeros((2, 2))
        bb[:, 0] = self.coordinates.min(1)
        bb[:, 1] = self.coordinates.max(1)
        return bb

    def area(self):
        raise NotImplementedError

    def length(self):
        raise NotImplementedError

    def sub_boundingbox(self, xedge=None, yedge=None):
        """
        Extract the x or y boundaries of the chaincode from a defined limits
        ``xedge`` or ``yedge``.

        Parameters
        ----------
        xedge : list, optional
            A list of two values to check.
        yedge : list, optional
            A list of two values to check.
        """
        if xedge is not None:
            edge = xedge
            IndexMask = 0
            IndexValue = 1
        elif yedge is not None:
            edge = yedge
            IndexMask = 1
            IndexValue = 0
        else:
            raise ValueError("Please input either `xedge` or `yedge`")
        mask = (self.coordinates[IndexMask, :] >= edge[0]) & (self.coordinates[IndexMask, :] <= edge[1])
        mx = np.ma.masked_array(self.coordinates[IndexValue, :], mask=(~mask))
        return [mx.min(), mx.max()]
