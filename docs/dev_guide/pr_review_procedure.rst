.. _pr_review:

**************************************
Pull Request Review and GitHub Members
**************************************

This document describes the procedure that the people in the "sunpy-maintainers" group must follow when merging a pull request.
If in any doubt about merging a pull request, the correct course of action is to contact the Lead Developer and ask.

Each pull request **must** meet the following criteria before it is considered for merge:

* The code must be PEP 8 compliant and meet the ill-defined SunPy quality standards.
  We have these in the :ref:`coding-standards` page.

* The PR must contain a CHANGELOG entry if it changes the behavior of any code.

* The test coverage should not decrease, and for new features should be at or very close to 100%.

* All code must be properly documented.
  Each function and each class must have an associated documentation string in the correct format.

Review Process
==============

Before the "merge" button is clicked the following criteria must be met:

* All the continuous integration must pass unless there is a known issue.

* At least two members (not the author of the PR) of the "sunpy-developers" group must have approved the PR, one should be a relevant subpackage maintainer.

* All comments posted on the thread must be resolved.

It is important that approval for merging the PR is done on the comment thread, as this becomes part of the "permanent record", this includes in during community meetings or in chat.

Continuous Integration
======================

Currently we have a variety of services that respond or activate on an opened pull request:

Comments from bots:

* `pep8speaks <https://github.com/OrkoHunter/pep8speaks>`_: Performs a PEP8 check on any submitted code. This is updated as the code changes.

Checks (they are located at the bottom of the pull request):

* `figure-tests (CircleCi) <https://circleci.com/gh/sunpy/sunpy/>`_: Runs two figure tests environments ("py38-figure", "py37-figure-devdeps").

* figure_report/figure_tests (Giles): Show the final results of the figure tests.

* figure_report_devdeps (Giles): Show the final results of the figure tests using development packages.

* changelog: absent | found (Giles): If a changelog is needed, this will check and will pass if a changelog with the correct number is found.

* milestone: absent | present (Giles): Will check that the pull request has a milestone assigned.

* `docs/readthedocs.org:sunpy (Read the Docs) <https://readthedocs.org/projects/sunpy/>`_: This builds our documentation.

* `sunpy.sunpy (Azure Pipelines) <https://dev.azure.com/sunpy/sunpy/_build>`_: Runs our test suite on all three operating systems.
  There are 10 separate checks for this.

* `codecov/patch (CodeCov) <https://codecov.io/gh/sunpy/sunpy/>`_: Checks how many lines of the code lack test coverage for the submitted code.

* `codecov/project (CodeCov) <https://codecov.io/gh/sunpy/sunpy/>`_: Checks how many lines of the code lack test coverage in sunpy overall.

It is common to see some of these checks fail.
This can be happen to due a change that has broken a test or a remote server that we use has failed.
Therefore it is important to check why a task failed and if has a pre-existing issue, it can be safe to ignore on that pull request.
However, you should try to ensure that as many checks pass before merging.

SunPy GitHub Groups
===================

This document has already referred to two SunPy groups, namely "developers" and "maintainers" there is also a third primary SunPy group "owners".
These owners' have control over the means of production.

SunPy owners
------------

The SunPy owners group is the group of people who have total control over the SunPy GitHub organization.
The SunPy board have control over who is in this group, it has been decided that generally it will be the Lead Developer and the SunPy board chair and vice-chair.

SunPy Maintainers
-----------------

This is the group of people who have push access to the main SunPy repository.
The membership of this group is at the discretion of the Lead Developer, but shall generally be made up of people who have demonstrated themselves to be trust worthy and active contributors to the project.

This group has `subgroups <https://github.com/orgs/sunpy/teams/sunpy-maintainers/teams>`__ for each section of the repository that has `maintainers <https://sunpy.org/team#maintainer-list>`__.
The members of these groups will automatically be requested to review all PRs which change files in that subpackage.



SunPy Developers
----------------

The members of this group have "read" access to the SunPy repository.
As all these repository are open anyway, what this effectively means is that these people can be assigned to issues.
The members of this group are people who are involved in the development of SunPy at a good frequency, they are people who’s opinions have been demonstrated to be constructive and informative.
