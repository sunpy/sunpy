"""
Modification of Chris Beaumont's mpl-modest-image package to allow the use of
set_extent.
"""
# This file is copied from glue under the terms of the 3 Clause BSD licence. See licenses/GLUE.rst

from __future__ import print_function, division

import matplotlib
rcParams = matplotlib.rcParams

import matplotlib.image as mi
import matplotlib.colors as mcolors
import matplotlib.cbook as cbook
from matplotlib.transforms import IdentityTransform, Affine2D

import numpy as np

IDENTITY_TRANSFORM = IdentityTransform()


class ModestImage(mi.AxesImage):

    """
    Computationally modest image class.

    ModestImage is an extension of the Matplotlib AxesImage class
    better suited for the interactive display of larger images. Before
    drawing, ModestImage resamples the data array based on the screen
    resolution and view window. This has very little affect on the
    appearance of the image, but can substantially cut down on
    computation since calculations of unresolved or clipped pixels
    are skipped.

    The interface of ModestImage is the same as AxesImage. However, it
    does not currently support setting the 'extent' property. There
    may also be weird coordinate warping operations for images that
    I'm not aware of. Don't expect those to work either.
    """

    def __init__(self, *args, **kwargs):
        self._pressed = False
        self._full_res = None
        self._full_extent = kwargs.get('extent', None)
        super(ModestImage, self).__init__(*args, **kwargs)
        self.invalidate_cache()
        self.axes.figure.canvas.mpl_connect('button_press_event', self._press)
        self.axes.figure.canvas.mpl_connect('button_release_event', self._release)
        self.axes.figure.canvas.mpl_connect('resize_event', self._resize)

        self._timer = self.axes.figure.canvas.new_timer(interval=500)
        self._timer.single_shot = True
        self._timer.add_callback(self._resize_paused)

    def remove(self):
        super(ModestImage, self).remove()
        self._timer.stop()
        self._timer = None

    def _resize(self, *args):
        self._pressed = True
        self._timer.start()

    def _resize_paused(self, *args):
        # If the artist has been removed, self.axes is no longer defined, so
        # we can return early here.
        if self.axes is None:
            return
        self._pressed = False
        self.axes.figure.canvas.draw_idle()

    def _press(self, *args):
        self._pressed = True

    def _release(self, *args):
        self._pressed = False
        self.stale = True
        self.axes.figure.canvas.draw_idle()

    def set_data(self, A):
        """
        Set the image array

        ACCEPTS: numpy/PIL Image A
        """
        self._full_res = A
        self._A = A

        if self._A.dtype != np.uint8 and not np.can_cast(self._A.dtype,
                                                         float):
            raise TypeError("Image data can not convert to float")

        if (self._A.ndim not in (2, 3) or
                (self._A.ndim == 3 and self._A.shape[-1] not in (3, 4))):
            raise TypeError("Invalid dimensions for image data")

        self.invalidate_cache()

    def invalidate_cache(self):
        self._bounds = None
        self._imcache = None
        self._rgbacache = None
        self._oldxslice = None
        self._oldyslice = None
        self._sx, self._sy = None, None
        self._pixel2world_cache = None
        self._world2pixel_cache = None

    def get_cursor_data(self, event):
        return None

    def contains(self, mouseevent):
        if self._A is None or self._A.shape is None:
            return False
        else:
            return super(ModestImage, self).contains(mouseevent)

    def set_extent(self, extent):
        self._full_extent = extent
        self.invalidate_cache()
        mi.AxesImage.set_extent(self, extent)

    def get_array(self):
        """Override to return the full-resolution array"""
        return self._full_res

    @property
    def _pixel2world(self):

        if self._pixel2world_cache is None:

            # Pre-compute affine transforms to convert between the 'world'
            # coordinates of the axes (what is shown by the axis labels) to
            # 'pixel' coordinates in the underlying array.

            extent = self._full_extent

            if extent is None:

                self._pixel2world_cache = IDENTITY_TRANSFORM

            else:

                self._pixel2world_cache = Affine2D()

                self._pixel2world.translate(+0.5, +0.5)

                self._pixel2world.scale((extent[1] - extent[0]) / self._full_res.shape[1],
                                        (extent[3] - extent[2]) / self._full_res.shape[0])

                self._pixel2world.translate(extent[0], extent[2])

            self._world2pixel_cache = None

        return self._pixel2world_cache

    @property
    def _world2pixel(self):
        if self._world2pixel_cache is None:
            self._world2pixel_cache = self._pixel2world.inverted()
        return self._world2pixel_cache

    def _scale_to_res(self):
        """
        Change self._A and _extent to render an image whose resolution is
        matched to the eventual rendering.
        """

        # Find out how we need to slice the array to make sure we match the
        # resolution of the display. We pass self._world2pixel which matters
        # for cases where the extent has been set.
        x0, x1, sx, y0, y1, sy = extract_matched_slices(axes=self.axes,
                                                        shape=self._full_res.shape,
                                                        transform=self._world2pixel)

        # Check whether we've already calculated what we need, and if so just
        # return without doing anything further.
        if (self._bounds is not None and
                sx >= self._sx and sy >= self._sy and
                x0 >= self._bounds[0] and x1 <= self._bounds[1] and
                y0 >= self._bounds[2] and y1 <= self._bounds[3]):
            return

        # Slice the array using the slices determined previously to optimally
        # match the display
        self._A = self._full_res[y0:y1:sy, x0:x1:sx]
        self._A = cbook.safe_masked_invalid(self._A)

        # We now determine the extent of the subset of the image, by determining
        # it first in pixel space, and converting it to the 'world' coordinates.

        # See https://github.com/matplotlib/matplotlib/issues/8693 for a
        # demonstration of why origin='upper' and extent=None needs to be
        # special-cased.

        if self.origin == 'upper' and self._full_extent is None:
            xmin, xmax, ymin, ymax = x0 - .5, x1 - .5, y1 - .5, y0 - .5
        else:
            xmin, xmax, ymin, ymax = x0 - .5, x1 - .5, y0 - .5, y1 - .5

        xmin, ymin, xmax, ymax = self._pixel2world.transform([(xmin, ymin), (xmax, ymax)]).ravel()

        mi.AxesImage.set_extent(self, [xmin, xmax, ymin, ymax])
        # self.set_extent([xmin, xmax, ymin, ymax])

        # Finally, we cache the current settings to avoid re-computing similar
        # arrays in future.
        self._sx = sx
        self._sy = sy
        self._bounds = (x0, x1, y0, y1)

        self.changed()

    def draw(self, renderer, *args, **kwargs):
        if self._full_res.shape is None:
            return
        if not self._pressed or self._bounds is None:
            self._scale_to_res()
        # Due to a bug in Matplotlib, we need to return here if all values
        # in the array are masked.
        if hasattr(self._A, 'mask') and np.all(self._A.mask):
            return
        super(ModestImage, self).draw(renderer, *args, **kwargs)


def main():
    from time import time
    import matplotlib.pyplot as plt
    x, y = np.mgrid[0:2000, 0:2000]
    data = np.sin(x / 10.) * np.cos(y / 30.)

    f = plt.figure()
    ax = f.add_subplot(111)

    # try switching between
    artist = ModestImage(ax, data=data)

    ax.set_aspect('equal')
    artist.norm.vmin = -1
    artist.norm.vmax = 1

    ax.add_artist(artist)

    t0 = time()
    plt.gcf().canvas.draw_idle()
    t1 = time()

    print("Draw time for %s: %0.1f ms" % (artist.__class__.__name__,
                                          (t1 - t0) * 1000))

    plt.show()


def imshow(axes, X, cmap=None, norm=None, aspect=None,
           interpolation=None, alpha=None, vmin=None, vmax=None,
           origin=None, extent=None, shape=None, filternorm=1,
           filterrad=4.0, imlim=None, resample=None, url=None, **kwargs):
    """Similar to matplotlib's imshow command, but produces a ModestImage

    Unlike matplotlib version, must explicitly specify axes
    """
    if norm is not None:
        assert(isinstance(norm, mcolors.Normalize))
    if aspect is None:
        aspect = rcParams['image.aspect']
    axes.set_aspect(aspect)
    im = ModestImage(axes, cmap=cmap, norm=norm, interpolation=interpolation,
                     origin=origin, extent=extent, filternorm=filternorm,
                     filterrad=filterrad, resample=resample, **kwargs)

    im.set_data(X)
    im.set_alpha(alpha)
    axes._set_artist_props(im)

    if im.get_clip_path() is None:
        # image does not already have clipping set, clip to axes patch
        im.set_clip_path(axes.patch)

    # if norm is None and shape is None:
    #    im.set_clim(vmin, vmax)
    if vmin is not None or vmax is not None:
        im.set_clim(vmin, vmax)
    # elif norm is None:
    #     im.autoscale_None()

    im.set_url(url)

    # update ax.dataLim, and, if autoscaling, set viewLim
    # to tightly fit the image, regardless of dataLim.
    if extent is not None:
        im.set_extent(extent)

    axes.add_image(im)

    def remove(h):
        axes.images.remove(h)

    im._remove_method = remove

    return im


def extract_matched_slices(axes=None, shape=None, extent=None,
                           transform=IDENTITY_TRANSFORM):
    """Determine the slice parameters to use, matched to the screen.

    :param ax: Axes object to query. It's extent and pixel size
               determine the slice parameters

    :param shape: Tuple of the full image shape to slice into. Upper
               boundaries for slices will be cropped to fit within
               this shape.

    :rtype: tulpe of x0, x1, sx, y0, y1, sy

    Indexing the full resolution array as array[y0:y1:sy, x0:x1:sx] returns
    a view well-matched to the axes' resolution and extent
    """

    # Find extent in display pixels (this gives the resolution we need
    # to sample the array to)
    ext = (axes.transAxes.transform([(1, 1)]) - axes.transAxes.transform([(0, 0)]))[0]

    # Find the extent of the axes in 'world' coordinates
    xlim, ylim = axes.get_xlim(), axes.get_ylim()

    # Transform the limits to pixel coordinates
    ind0 = transform.transform([min(xlim), min(ylim)])
    ind1 = transform.transform([max(xlim), max(ylim)])

    def _clip(val, lo, hi):
        return int(max(min(val, hi), lo))

    # Determine the range of pixels to extract from the array, including a 5
    # pixel margin all around. We ensure that the shape of the resulting array
    # will always be at least (1, 1) even if there is really no overlap, to
    # avoid issues.
    y0 = _clip(ind0[1] - 5, 0, shape[0] - 1)
    y1 = _clip(ind1[1] + 5, 1, shape[0])
    x0 = _clip(ind0[0] - 5, 0, shape[1] - 1)
    x1 = _clip(ind1[0] + 5, 1, shape[1])

    # Determine the strides that can be used when extracting the array
    sy = int(max(1, min((y1 - y0) / 5., np.ceil(abs((ind1[1] - ind0[1]) / ext[1])))))
    sx = int(max(1, min((x1 - x0) / 5., np.ceil(abs((ind1[0] - ind0[0]) / ext[0])))))

    return x0, x1, sx, y0, y1, sy


if __name__ == "__main__":
    main()
