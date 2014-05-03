"""
DOES STUFF TO YOUR IMAGE.

IT WILL NEVER BE THE SAME EVERY AGAIN (UNLESS YOU MAKE A COPY)
"""

from __future__ import division
import itertools
import numpy as np
import scipy.ndimage as ndimage


def multiscale_gaussian(data, sigma=[1.25, 2.5, 5, 10, 20, 40], k=0.7,
                        gamma=3.2, h=0.7, weights=None):
    """
    Multi-scale Gaussian normalization.

    Parameters
    ----------
    data : numpy.ndarray
        Image to be transformed.

    sigma : list, optional
        Range of guassian widths to transform over.

    k : float, optional
        Controls the severity of the arctan transformation.

    gamma : float, optional
        The value used to calulcate the  global gamma-transformed image.
        Ideally should be betwen 2.5 to 4.

    h : float, optional
        Global weight.

    weights : list, optional
        Used to weight all the transformed images during the calculation of the
        final image. If not specificed, all weights are one.

    Returns
    -------
    image: numpy.ndarray
        Normalized image.

    Reference
    ---------
    Morgan, Huw, and Miloslav Druckmuller.
    "Multi-scale Gaussian normalization for solar image processing."
    arXiv preprint arXiv:1403.6613 (2014).

    Notes
    -----
    In practice, the weights and h may be adjusted according to the desired
    output, and also according to the type of input image
    (e.g. wavelength or channel).
    For most purposes, the weights can be set
    equal for all scales.
    """
    if not weights:
        weights = np.ones(len(sigma))

    data[data <= 0] = 1e-15  # Makes sure that all values are above zero
    image = np.zeros(data.shape)

    for s, weight in itertools.izip(sigma, weights):
        conv = ndimage.filters.gaussian_filter(data, sigma=s)
        lm_sub = data - conv
        sigmaw = np.sqrt(ndimage.filters.gaussian_filter(lm_sub ** 2, sigma=s))
        image += np.arctan(k * (lm_sub / sigmaw)) * weight

    data_min = data.min()
    data_max = data.max()
    global_gamma = h * ((data-data_min)/(data_max-data_min)) ** (1/gamma)

    return global_gamma + ((1-h)/len(sigma)) * image