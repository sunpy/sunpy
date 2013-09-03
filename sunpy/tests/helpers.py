# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 18:04:34 2013

@author: stuart
"""
import tempfile

import numpy as np
import matplotlib.pyplot as plt

def plot_comparison(fig, filename):
    tmp = tempfile.NamedTemporaryFile()
    fig.savefig(tmp, format='png')
    fname = tmp.name
    tmp.delete = False
    tmp.close()
    im1 = plt.imread(filename)
    im2 = plt.imread(fname, format='png')
    pass_ = np.allclose(im1, im2, rtol=1e-3)
    if pass_:
        return pass_
    else:
        raise Exception("Plots do not match")
