from __future__ import division

import numpy as np
import scipy.ndimage as ndimage


def multiscale_gaussian(data, sigma=[1.25, 2.5, 5, 10, 20, 40], k=0.7,
                        gamma=3.2, h=0.7, weights=None, truncate=3):
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
        Ideally should be between 2.5 to 4.

    h : float, optional
        Global weight.

    weights : list, optional
        Used to weight all the transformed images during the calculation of the
        final image. If not specificed, all weights are one.

    width : `int`
        An odd integer defining the width of the kernel to be convolved.

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

    # 1. Replace spurious negative pixels with zero
    data[data <= 0] = 1e-15  # Makes sure that all values are above zero
    image = np.empty(data.shape, dtype=data.dtype)
    conv = np.empty(data.shape, dtype=data.dtype)
    sigmaw = np.empty(data.shape, dtype=data.dtype)

    for s, weight in zip(sigma, weights):
        # 2 & 3 Create kernel and convolve with image
        ndimage.filters.gaussian_filter(data, sigma=s,
                                        truncate=truncate, output=conv)
        # 5. Calculate difference between image and the local mean image,
        # square the difference, and convolve with kernel. Square-root the
        # resulting image to give ‘local standard deviation’ image sigmaw
        conv -= data
        conv **= 2
        ndimage.filters.gaussian_filter(conv, sigma=s,
                                        truncate=truncate, output=sigmaw)
        np.sqrt(sigmaw, out=sigmaw)
        conv /= sigmaw

        # 6. Apply arctan transformation on Ci to give C'i
        conv *= k
        conv *= weight
        np.arctan(conv, out=conv)

        image += conv

    # delete these arrays here as it reduces the total memory consumption when
    # we create the Cprime_g temp array below.
    del conv
    del sigmaw

    # 8. Take weighted mean of C'i to give a weighted mean locally normalised
    # image.
    image /= len(sigma)

    # 9. Calculate global gamma-transformed image C'g
    data_min = data.min()
    data_max = data.max()
    Cprime_g = (data - data_min)
    Cprime_g /= (data_max - data_min)
    Cprime_g **= (1/gamma)
    Cprime_g *= h

    image *= (1 - h)
    image += Cprime_g

    return image
