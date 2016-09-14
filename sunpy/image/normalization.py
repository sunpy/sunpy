from __future__ import division
import itertools
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
    image = np.zeros(data.shape)

    for s, weight in zip(sigma, weights):
        # 2 & 3 Create kernel and convolve with image
        conv = ndimage.filters.gaussian_filter(data, sigma=s, truncate=truncate)
        # 5. Calculate difference between image and the local mean image,
        # square the difference, and convolve with kernel. Square-root the
        # resulting image to give ‘local standard deviation’ image sigmaw
        lm_sub = data - conv
        np.subtract(data, conv, out=conv)
        lm_sub = conv
        sigmaw = np.sqrt(ndimage.filters.gaussian_filter(lm_sub ** 2, sigma=s, truncate=truncate))
        Ci = lm_sub / sigmaw
        # 6. Apply arctan transformation on Ci to give C'i
        image += np.arctan(k * Ci) * weight

    # 8. Take weighted mean of C'i to give a weighted mean locally normalised
    # image.
    image /= len(sigma)

    # 9. Calculate global gamma-transformed image C'g
    data_min = data.min()
    data_max = data.max()
    Cprime_g = ((data - data_min) / (data_max - data_min)) ** (1/gamma)

    I = (h * Cprime_g) + ((1 - h) * image)

    return I
