from __future__ import absolute_import
from __future__ import division

import sys

__all__ = ['TTYProgressBar']

class TTYProgressBar(object):
    """
    A simple progress bar to visualize progress on a TTY (teletypewriter).
    It needs to support '\b' to delete characters.

    The ProgressBar interface is start, finish, draw and poke.
    """
    SYMBOL = '='
    LEFT_BORDER = '['
    RIGHT_BORDER = ']'

    def __init__(self, n, current=0, width=40, output=sys.stdout):
        """
        Parameters
        ----------
        n : int
            Total number of items until completion
        current : int
            Current state of completion
        width : int
            Width of progress bar.
        output : file
            Teletypewriter to print on.
        """
        self.n = n
        self.current = current
        self.width = width
        # Future division, so it is float.
        self.step = self.n / self.width
        self.output = output

    def start(self):
        """
        Draw empty bar to output.
        """
        self.output.write(
            self.LEFT_BORDER + " " * (len(self.SYMBOL) * self.width) +
                self.RIGHT_BORDER
        )
        self.output.flush()
        self.output.write("\b" * (self.width+len(self.RIGHT_BORDER)))

    def finish(self):
        """
        Finish the bar, the ProgressBar cannot be used after this
        method was called.
        """
        print

    def _draw_one(self):
        """
        Advance progress bar by one.
        """
        self.output.write(self.SYMBOL)

    def draw(self):
        """
        Draw current state of progress bar onto and empty bar.
        """
        cur = self.current
        self.current = 0
        for _ in xrange(cur):
            self.poke()

    def poke(self, n=1):
        """
        Increase finished items by n. May advance the progress bar by one
        or more fields.
        """
        if self.current > self.n:
            raise ValueError("ProgressBar overflowed.")

        diff = int((self.current + n) / self.step) - int(self.current / self.step)
        for _ in xrange(diff):
            self._draw_one()
        self.current += n
