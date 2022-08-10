.. _pr_checklist:

***********************
Pull Request Check List
***********************

The pull request (commonly referred to as a PR) check list below is an outline of the steps that should be taken when making a contribution to a SunPy repository on Github.

New to Git and Github?
Check out the :ref:`newcomers` and `Git Cheat Sheets <https://training.github.com/>`__.

If you would prefer a visual Git interface, you can try `Github Desktop <https://desktop.github.com/>`__ or `GitKraken <https://www.gitkraken.com/>`__.

#. Review and test changes locally on your machine (see :ref:`testing`).
    #. Double check that a pull request does not exist for the changes you are making.
       Ideally, check that there is an issue that details what you want to change and why.
    #. If you are contributing code, review the :ref:`coding-standards` page.
    #. See the :ref:`dev_guide` for guidelines regarding code tests, documentation, and other types of contributions.
    #. Have you tested your changes with the latest version of ``sunpy``?
       If not, update your local copy from your `remote repository <https://docs.github.com/en/get-started/using-git/getting-changes-from-a-remote-repository>`__ on Github.
#. Push your changes to Github.
    #. `Create a new branch <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-and-deleting-branches-within-your-repository>`__ in your fork of ``sunpy``.
    #. Give this branch a name that reflects the changes you are making.
    #. Create commits that describe the changes.
    #. `Push your changes <https://docs.github.com/en/get-started/using-git/pushing-commits-to-a-remote-repository>`__ to the new branch on your remote fork.
#. `Compare your branch <https://docs.github.com/en/pull-requests/committing-changes-to-your-project/viewing-and-comparing-commits/comparing-commits>`__ with ``sunpy/main``.
    #. Resolve any conflicts that occur ideally with a `git rebase <https://www.atlassian.com/git/tutorials/rewriting-history/git-rebase>`__.
#. Create a pull request.
    #. `Create a pull request <https://docs.github.com/en/get-started/quickstart/hello-world#opening-a-pull-request>`__ from the branch you have created/updated.
    #. Add a description to the pull request that describes the changes you have made.
       Remember to delete the preamble within the message box.
    #. `Link <https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/autolinked-references-and-urls>`__ to any relevant issues, pull requests, or comments in your description.
#. Add a changelog to your pull request.
    #. A `changelog <https://github.com/sunpy/sunpy/tree/main/changelog#changelog>`__ is a short record of the type of changes made in your pull request.
       Other users are the intended audience, and you can have multiple logs per pull request.
#. Maintainers will review your pull request :ref:`pr_review`.
    #. Tweak anything that others highlight and push the changes to your branch.
       You can also commit suggestions either in bulk or single commits via the Github user interface.
    #. Discuss possible changes or improvements in the comments with the reviewers.
    #. Review the Continuous Integration (CI) :ref:`ci_jobs` tests and fix any errors or warnings that are found.
        #. If you are confused by an error that the continuous integration is giving you, submit a comment in your pull request.
#. Ask questions if you get stuck or confused at any point!
    #. Open-source projects are about communication and collaboration.
    #. `Join the SunPy Matrix Chat channel <https://matrix.to/#/+sunpy:openastronomy.org>`__.

This guide is partially based on Astropy's `Development Workflow <https://docs.astropy.org/en/latest/development/workflow/development_workflow.html>`__.
