# -*- coding: utf-8 -*-

import warnings

from sunpy.util.decorators import deprecated
from sunpy.util.decorators import add_common_docstring
from sunpy.visualization.animator.mapsequenceanimator import MapSequenceAnimator

__all__ = ['MapCubeAnimator']


@deprecated('0.9.1', message='MapCubeAnimator deprecated in favor of MapSequenceAnimator. MapSequence has the same functionality as MapCube.')
@add_common_docstring(append='MapCubeAnimator deprecated in favor of MapSequenceAnimator. MapSequence has the same functionality as MapCube.',
                      prepend=MapSequenceAnimator.__doc__)
def MapCubeAnimator(*args, **kwargs):
    return MapSequenceAnimator(*args, **kwargs)
