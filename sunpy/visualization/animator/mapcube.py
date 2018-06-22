# -*- coding: utf-8 -*-

from sunpy.util.decorators import deprecated
from sunpy.visualization.animator.mapsequenceanimator import MapSequenceAnimator

__all__ = ['MapCubeAnimator']


MapCubeAnimator = deprecated(MapSequenceAnimator, "0.9.1",
                             message='MapCubeAnimator deprecated in favor of MapSequenceAnimator. \
                                      MapSequence has the same functionality as MapCube.')
