# -*- coding: utf-8 -*-

import warnings

from sunpy.util.exceptions import SunpyDeprecationWarning
from sunpy.visualization.mapsequenceanimator import MapSequenceAnimator

__all__ = ['MapCubeAnimator']

warnings.warn("Deprecated in favor of MapSequenceAnimator. MapSequence has the same functionality as MapCube.",
              SunpyDeprecationWarning, stacklevel=2)

MapCubeAnimator = MapSequenceAnimator
