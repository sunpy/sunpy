import sunpy.visualization.animator
from sunpy.util.decorators import deprecated

__all__ = ['MapSequenceAnimator']


@deprecated("5.0", alternative="sunpy.visualization.animator.MapSequenceAnimator")
class MapSequenceAnimator(sunpy.visualization.animator.MapSequenceAnimator):
    pass
