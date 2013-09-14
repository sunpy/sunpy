# -*- coding: utf-8 -*-
# Author:   Michael Malocha <mjm159@humboldt.edu>
# Last Edit:  September 14th, 2013
#
# This module was developed with funding from the GSOC 2013 summer of code
#

from __future__ import absolute_import
import sys

__author__ = 'Michael Malocha'
__version__ = 'September 14th, 2013'


def progress_bar(phrase, position, total, bar_size=20):
    """
    Prints a simple progress bar to the screen

    Parameters
    ----------
    phrase: str
        The desired phrase to precede the progress bar with each update.
    position: int
        The current location in the data set (will be divided by 'total'
        to determine percent finished).
    total: int
        The total size of the data set
    bar_size: int
        Size of the bar, defaults to 20 characters in length

    Examples
    --------
    >>> progress_bar("Progress:", 10, 100)
    Progress: [##                  ] 10%

    >>> progress_bar("Progress:", 35, 100)
    Progress: [#######             ] 35%

    >>> progress_bar("Progress:", 35, 83)
    Progress: [########            ] 42%
    """
    position = float(format(position, '.1f'))
    fraction = position / total
    percent = str(int(100 * fraction)) + '%'
    place = int(fraction * bar_size)
    pounds = '#' * place
    blanks = ' ' * (20 - place)
    prog_bar = ' [' + pounds + blanks + '] ' + percent
    sys.stdout.write('\r' + ' ' * 52)
    sys.stdout.flush()
    sys.stdout.write('\r' + phrase + prog_bar)
    sys.stdout.flush()
