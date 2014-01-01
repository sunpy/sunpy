# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.animation as anim
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

def toggle_pylab(fn):
    """ A decorator to prevent functions from opening matplotlib windows
        unexpectedly when sunpy is run in interactive shells like ipython
        --pylab.

        Toggles the value of matplotlib.pyplot.isinteractive() to preserve the
        users' expections of pylab's behaviour in general. """

    if plt.isinteractive():
        def fn_itoggle(*args, **kwargs):
            plt.ioff()
            ret = fn(*args, **kwargs)
            plt.ion()
            return ret
        return fn_itoggle
    else:
        return fn

def animate_array(array, n, axes=None, cmap=plt.get_cmap('gray'), norm='static',
                  extent=None, interval=200, colorbar=False, ani_args={},
                  **imshow_args):
    """
    Create a mpl animation along the ith axis of an array

    Parameters
    ----------

    norm: 'dynamic', 'static', or mpl.Normalize, or list of mpl.Normalize
    """

    def get_slice(i, n, nax=3):
        """
        Slice an array, i'th element on n'th axes

        Parameters
        ----------
        i: int
            The element to select
        n: int
            The axis along which to index the i'th element
        nax: int
            The total number of axes in the array
        """
        arr_slice = [slice(None)]*nax
        arr_slice[n] = i
        return arr_slice

    n_im = array.shape[n]

    if not axes:
        axes = plt.gca()
    fig = axes.get_figure()

    if isinstance(norm, basestring):
        if norm == 'static':
            norm = [Normalize(vmin=array[get_slice(0,n)].min(),
                              vmax=array[get_slice(0,n)].max())]*n_im
        elif norm == 'dynamic':
            norm = []
            for j in xrange(0,n_im):
                norm.append(Normalize(vmin=array[get_slice(j,n)].min(),
                                      vmax=array[get_slice(j,n)].max()))
        else:
            raise TypeError("norm is not my cow")

    elif isinstance(norm, Normalize):
        norm = [norm]*n_im

    elif not isinstance(norm, list):
        raise TypeError("norm is not my cow")

    #make imshow kwargs a dict
    kwargs = {'origin':'lower',
              'cmap':cmap,
              'norm':norm[0],
              'extent':extent,
              'interpolation':'nearest'}
    kwargs.update(imshow_args)

    im = axes.imshow(array[get_slice(0,n)], **kwargs)

    #Set current image (makes colorbar work)
    plt.sci(im)

    if colorbar:
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="5%", pad=0.2)
        cbar = plt.colorbar(im,cax)

    def updatefig(i, *args):
        """
        args is [im,array, n]
        """
        im = args[0]
        im.set_array(args[2][get_slice(i,n)])
        im.set_norm(norm[i])
#        if args[1]:
#            axes.set_title("%s %s" % (self[i].name, self[i].date))
        axes.set_title("%i"%i)
    ani = anim.FuncAnimation(fig, updatefig,
                                        frames=xrange(0,n_im),
                                        fargs=[im,False,array,n],
                                        interval=interval,
                                        blit=False, **ani_args)

    return ani