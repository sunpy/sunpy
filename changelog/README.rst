.. This README was adapted from the pytest changelog readme under the terms of the MIT licence.

This directory contains "newsfragments" which are short files that contain a small **ReST**-formatted
text that will be added to the next ``CHANGELOG``.

The ``CHANGELOG`` will be read by users, so this description should be aimed to sunpy users
instead of describing internal changes which are only relevant to the developers.

Make sure to use full sentences with correct case and punctuation, for example::

    Add support for Helioprojective coordinates in `sunpy.coordinates.frames`.

Each file should be named like ``<PULL REQUEST>.<TYPE>.rst``, where
``<PULL REQUEST>`` is a pull request number, and ``<TYPE>`` is one of:

* ``breaking``: A change which requires users to change code and is not backwards compatible. (Not to be used for removal of deprecated features.)
* ``feature``: new user facing features and new behavior.
* ``bugfix``: fixes a reported bug.
* ``doc``: documentation improvement, like rewording an entire session or adding missing docs.
* ``removal``: feature deprecation or removal.
* ``trivial``: A change which has no user facing effect or is tiny.

So for example: ``123.feature.rst``, ``456.bugfix.rst``.

If you are not sure what pull request type to use, don't hesitate to ask in your PR.

Note that the ``towncrier`` tool will automatically
reflow your text, so it will work best if you stick to a single paragraph, but multiple sentences and links are OK
and encouraged. You can install ``towncrier`` and then run ``towncrier --draft``
if you want to get a preview of how your change will look in the final release notes.
