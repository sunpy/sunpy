.. _backports:

***********************************
Making Changes to Released Versions
***********************************

When changes need to be made in a bugfix release of an already released version of sunpy this is done by a process called "backporting".
The process is as follows:

* Open a Pull Request (PR) as normal targeting the ``main`` branch.
* Apply a label to that PR (or get a maintainer to do it) of the format ``Backport X.Y``.
* Once the PR has been merged a bot will open a new PR with the same changes targeting the ``X.Y`` branch.

Managing backports is done by the sunpy maintainers, as a contributor you generally shouldn't have to worry about it.
The following documentation is for maintainers on how to interact with the backport bot.


Controlling the Backport Bot
============================

The backport bot in use on the sunpy repo is `MeeseeksDev <https://github.com/MeeseeksBox/MeeseeksDev/>`__.

Upon merging a PR the bot will look for the existence of a label with the following template in the label description field ``on-merge: backport to X.Y``.

If the decision to backport a PR is taken after the merge of the PR, then a command needs to be added to a comment on the PR: ``@MeeseeksDev backport [to] {branch}``


Manual Backports
================

If a backport fails, meeseeks will add a comment to the PR and add a label named "Still Needs Manual Backport".
If you then manually backport the PR, please remove this label when the backport PR is open so that the label is an accurate list of PRs in need of manual backporting.

When doing a manual backport please do not forget to put the PR number you are backporting **in the description** of the PR so that GitHub links the two PRs together.
(This automatic linking does not happen if the number is in the title).
