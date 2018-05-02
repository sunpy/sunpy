.. _maintainer-workflow:

Workflow for Maintainers
========================

This page is for maintainers those of us who merge our own or other peoples' changes into the upstream repository.

Being as how you're a maintainer, you are completely on top of the basic stuff in :ref:`version_control`.

Integrating changes
-------------------

Whenever possible, merge pull requests automatically via the pull request manager on GitHub.
Merging should only be done manually if there is a really good reason to do this!

Please try to make sure that pull requests do not contain a messy history with merges, etc.
If this is the case, then you have two options.
To use the rebase merge on Github or to follow the manual instructions, and make sure the fork is rebased to tidy the history before committing.
The second option should only be used on complex pull requests.

Using Milestones and Labels
===========================

These guidelines are adapted from `similar guidelines <http://docs.astropy.org/en/stable/development/workflow/maintainer_workflow.html#using-milestones-and-labels>`_ followed by Astropy:

* All open pull requests should have a milestone.

* Only confirmed issues that are release critical or for some other reason should be addressed for a release, should have a milestone.

* In general there should be the following open milestones:

  * The next bug fix releases for any still-supported version lines; for example if 0.4 is in development and
    0.2.x and 0.3.x are still supported there should be milestones for the next 0.2.x and 0.3.x releases.

  * The next X.Y release, i.e. the next minor release; this is generally the next release that all development in
    master is aimed toward.

  * The next X.Y release +1; for example if 0.3 is the next release, there should also be a milestone for 0.4 for
    issues that are important, but that we know won't be resolved in the next release.

* When in doubt about which milestone to use for an issue, use the next minor release, it can always be moved once
  it's been more closely reviewed prior to release.

* Issues that require fixing in the mainline, but that also are confirmed to apply to supported stable version lines
  should be marked with a ``Affects Release`` label and the corresponding supported stable version label ``v0.4.x``.

Using Projects
==============

Projects allow us to layout current pull requests and issues in a manner that enables a more `meta` view regarding major releases.
We categorize pull requests and issues into several levels of priorities and whether these can be classed as blockers before a release can be attempted.
Further we can add general notes that someone deems important for a release.

Updating and Maintaining the Changelog
======================================

The SunPy "changelog" is kept in the file ``CHANGELOG.rst`` at the root of the repository.
As the filename extension suggests this is a reStructured Text file.
The purpose of this file is to give a technical, but still user (and developer) oriented overview of what changes were made to Sunpy between each public release.
The idea is that it's a little more to the point and easier to follow than trying to read through full git log.
It lists all new features added between versions, so that a user can easily find out from reading the changelog when a feature was added.
Likewise it lists any features or APIs that were changed (and how they were changed) or removed.
It also lists all bug fixes.
Affiliated packages are encouraged to maintain a similar changelog.

Adding to the changelog
-----------------------

We want to support a single method that one may take to adding a new entry to the changelog.
It should be said that *all* additions to the changelog should be made first in the 'master' branch.
This is because every release of Sunpy includes a copy of the changelog, and it should list all the changes in every prior version of Sunpy.
For example, when Sunpy v0.3.0 is released, in addition to the changes new to that version the changelog should have all the changes from every v0.2.x version (and earlier) released up to that point.

We want to support one approach for including a changelog entry for a new feature or bug fix:

* Include the changelog update in the same pull request as the change.
  That is, assuming this change is being made in a pull request it can include an accurate changelog update along with it.

  Pro: An addition to the changelog is just like any other documentation update, and should be part of any atomic change to the software.
  It can be pulled into master along with the rest of the change.

  Con: If many pull requests also include changelog updates, they can quickly conflict with each other and require rebasing.
  This is not difficult to resolve if the only conflict is in the changelog, but it can still be trouble especially for new contributors.

As any SunPy maintainer can push to any or all pull requests, this can be done as the last commit to to a pull request with a `[ci skip]` message once the pull request is ready to merge.

Changelog format
----------------

The exact formatting of the changelog content is a bit loose for now (though it might become stricter if we want to develop more tools around the
changelog).
The format can be mostly inferred by looking at previous versions.
Each release gets its own heading (using the ``=`` heading marker) with the version and release date.
Releases still under development have ``(unreleased)`` as there is no release date yet.

There are two subheadings (using the ``-`` marker): "New Features" and "Bug Fixes".
In case of other changes, "Other Changes and Additions" can be added (but is not default), mostly a catch-all for miscellaneous changes, though there's no reason not to make up additional sub-headings if it seems appropriate.

The actual texts of the changelog entries are typically just one to three sentences, they should be easy to glance over.
Most entries end with a reference to an issue/pull request number in square brackets.

A single changelog entry may also reference multiple small changes.
For example::

  - Minor documentation fixes and restructuring. [#935, #967, #978, #1004, #1028, #1047]

Beyond that, the best advice for updating the changelog is just to look at existing entries for previous releases and copy the format.
