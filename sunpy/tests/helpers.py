# -*- coding: utf-8 -*-

import tempfile
import os.path
import shutil

import numpy as np
import matplotlib.pyplot as plt

import sunpy

def plot_comparison(fig, filename):
    """
    A Function for testing plotting

    This function will compare a figure instance to a saved png file.
    the png file to compare to must be saved into 
    ```sunpy/tests/saved_plots/```
    the path to this dir will be added to the filename supplied.    
    
    Parameters
    ----------
    fig: mpl.figure
        A Figure instance to compare to
    
    filename: string
        The filename of the plot, without path
    """
    rootdir = os.path.join(os.path.dirname(sunpy.__file__), "tests", "saved_plots")
    filename = os.path.join(rootdir, filename)
    
    tmp = tempfile.NamedTemporaryFile()
    fig.savefig(tmp, format='png')
    fname = tmp.name
    tmp.delete = False
    tmp.close()

    # If provided file exists then compare the plot, else save the file to 
    # the current directory
    if os.path.exists(filename):
        im1 = plt.imread(filename)
        im2 = plt.imread(fname, format='png')
        pass_ = np.allclose(im1, im2, rtol=1e-3)
        if pass_:
            return pass_
        else:
            raise Exception("Plots do not match")
    else:
        shutil.copyfile(fname, "./" + os.path.basename(filename))
        raise Exception("File did not exist, copied to working directory")