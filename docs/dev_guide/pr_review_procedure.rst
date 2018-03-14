.. _pr_review:

This document describes the procedure that the people in the sunpy-maintainters group must follow when merging a pull request.
If in any doubt about merging a pull request, the correct course of action is to contact the Lead Developer and ask.

Fundamental Components of a Pull Request
========================================

Each pull request *must* meet the following criteria before it is considered for merge:

*  The code must be PEP 8 compliant and meet the ill-defined SunPy quality standards.
   We have these in the :ref:`coding-standards` page.
*  The PR must contain a CHANGELOG entry if it changes the behavior of any code.
*  The test coverage should not decrease, and for new features should be at or very close to 100%.
*  All code must be properly documented.
   Each function and each class must have an associated documentation string in the correct format.

Review Process
==============

Before the ‘merge’ button is clicked the following criteria must be met:

*  All the continuous integration pass unless there is a known issue.
*  At least two members (not the author of the PR) of the sunpy-developers group have approved the PR.
*  All comments posted on the thread must be resolved.

It is important that approval for merging the PR is done on the comment thread, as this becomes part of the ‘permanent record’, this includes in during community meetings or in chat.

Continuous Integration
======================

Currently we have a variety of bots or services that respond or activate on an opened pull request.
While we try not to change them, they have undergone several changes with the aim of making them clear and focused on specific issues.

*  pep8speaks: Performs a PEP8 check on any submitted code.
*  CircleCi: Tests to see if sunpy installs and builds the documentation.
*  Giles: Returns a link if the documentation does build successfully.
*  Travis: Runs our test suite to make sure it passes on Linux and mac OS.
*  AppVeyor: Runs our test suite to make sure it passes on Windows.
*  CodeCov: Checks how many lines of the code lack test coverage.

We support several custom tags you can add anywhere in the commit message.
Please use these tags extensively, especially for documentation PRs and WIP commits.

*  [skip ci] or [ci skip]  - Skips Appveyor and Travis.
*  [skip appveyor] - Skips Appveyor only.
*  [skip travis] or [travis skip] -  These terminate the Travis builds before the install steps, effectively skipping testing.
*  [docs only] or [build docs] - These terminate all non documentation build jobs in Travis before the install steps.

We have auto-cancellation enabled both on Appveyor and Travis for SunPy core.
This means that queued builds for commits are cancelled if there is a newer commit pushed to that given branch.

SunPy GitHub Groups
===================

This document has already referred to two SunPy groups, namely ‘developers’ and ‘maintainers’ there is also a third primary SunPy group ‘owners’.
These owners' have control of the means of production.

SunPy owners
------------

The SunPy owners group is the group of people who have total control over the SunPy GitHub organization.
The SunPy board have control over who is in this group, it has been decided that generally it will be the Lead Developer and the SunPy board chair and vice-chair.

SunPy Maintainers
-----------------

This is the group of people who have push access to the main SunPy repos.
The membership of this group is at the discretion of the Lead Developer, but shall generally be made up of people who have demonstrated themselves to be trust worthy and active contributors to the project.

SunPy Developers
----------------

The members of this group have ‘read’ access to the SunPy orgs repos.
As all these repos are open anyway, what this effectively means is that these people can be assigned to issues.
The members of this group are people who are involved in the development of SunPy at a good frequency, they are people who’s opinions have been demonstrated to be constructive and informative.

.. _review: https://help.github.com/articles/about-pull-request-reviews/
