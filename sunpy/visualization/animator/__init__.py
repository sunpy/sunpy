import mpl_animators as _animators

from sunpy.util.decorators import deprecated
from sunpy.visualization.animator.mapsequenceanimator import MapSequenceAnimator

__all__ = ['LineAnimator', 'ImageAnimator', 'ArrayAnimatorWCS',
           'ArrayAnimator', 'BaseFuncAnimator', 'MapSequenceAnimator']

deprecated_animator = deprecated("3.1", alternative="the new mpl-animators package")


@deprecated_animator
class BaseFuncAnimator(_animators.BaseFuncAnimator):
    pass


@deprecated_animator
class ArrayAnimator(_animators.ArrayAnimator):
    pass


@deprecated_animator
class ArrayAnimatorWCS(_animators.ArrayAnimatorWCS):
    pass


@deprecated_animator
class ImageAnimator(_animators.ImageAnimator):
    pass


@deprecated_animator
class LineAnimator(_animators.LineAnimator):
    pass
