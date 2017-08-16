.. _reporting-problems:

Reporting Bugs
==============

If you are having a problem with sunpy, search the mailing
lists first: it is possible that someone else has already run into
your problem.

If not, please provide the following information in your e-mail to the
`mailing list <http://groups.google.com/forum/#!forum/sunpy>`_:

  * your operating system; (Linux/UNIX users: post the output of ``uname -a``)

  * sunpy version::

        >>> import sunpy   # doctest: +SKIP
        >>> sunpy.util.system_info()   # doctest: +SKIP

  * how you obtained sunpy.

  * any customizations to your ``sunpyrc`` file (see
    :ref:`customizing-sunpy`).

  * Please try to provide a *minimal*,
    standalone Python script that demonstrates the problem.  This is
    *the* critical step.  If you can't post a piece of code that we
    can run and reproduce your error, the chances of getting help are
    significantly diminished.  Very often, the mere act of trying to
    minimize your code to the smallest bit that produces the error
    will help you find a bug in *your* code that is causing the
    problem.

You will likely get a faster response writing to the mailing list than
filing a bug in the `bug tracker <http://github.com/sunpy/sunpy/issues>`_.
If your problem has been determined to be a bug and can not be quickly solved, the issues
may be filed a bug in the tracker so the issue doesn't get lost.
