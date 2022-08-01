.. _cheatsheet:

******************************
SunPy Pull Request Cheat Sheet
******************************

The pull request (PR) cheat sheet below is a rough outline of the steps that should be taken when making a contribution to SunPy on Github.

1. Review and test changes locally on your machine
    a. Double check that an Issue or PR does not exist for the changes you are making
    b. If you are contributing code, review the `Coding Standards page <https://docs.sunpy.org/en/latest/dev_guide/contents/code_standards.html>`_
    c. See the `Developer's Guide <https://docs.sunpy.org/en/latest/dev_guide/index.html>`_ for guidelines regarding code tests, documentation, and other types of contributions
    d. Have you tested your changes with the latest version of SunPy? If not, update your local copy of SunPy from your `remote repository <https://docs.github.com/en/get-started/using-git/getting-changes-from-a-remote-repository>`_ on Github
2. Push your changes to Github
    a. `Create a new branch <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-and-deleting-branches-within-your-repository>`_ of sunpy/main in your fork of SunPy
    b. Give this branch a name that reflects the changes you are making
    c. Make multiple commits that describe the changes
    d. `Push your changes <https://docs.github.com/en/get-started/using-git/pushing-commits-to-a-remote-repository>`_ to the new branch
3. `Compare your branch <https://docs.github.com/en/pull-requests/committing-changes-to-your-project/viewing-and-comparing-commits/comparing-commits>`_ with sunpy/main
    a. Resolve any conflicts that occur
4. Create a pull request
    a. `Create a pull request <https://docs.github.com/en/get-started/quickstart/hello-world#opening-a-pull-request>`_ from the branch you have updated
    b. Add a description to the PR that describes the changes you have made
    c. `Link <https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/autolinked-references-and-urls>`_ to any relevant Issues or PRs in your description
5. Add a changelog to your PR
    a. A `changelog <https://github.com/sunpy/sunpy/tree/main/changelog#changelog>`_ is a short record of the type of changes made in your PR
6. SunPy maintainers will `review your PR <https://docs.sunpy.org/en/latest/dev_guide/contents/pr_review_procedure.html>`_
    a. Fix any errors that others highlight and push the changes to your branch
    b. Discuss possible changes/improvements in the comments
    c. Review the `Continuous Integration (CI) <https://docs.sunpy.org/en/latest/dev_guide/contents/ci_jobs.html>`_ tests and fix any errors that are found
        i. If you are confused by an error that the CI is giving you just ask about in your PR
7. Ask questions if you get stuck or confused at any point
    a. Open-source projects are about communication and collaboration
    b. Join the SunPy Matrix Chat channel (link at the bottom of this page)

Making a change to the SunPy website? See `these guidelines <https://github.com/sunpy/sunpy.org#sunpyorg-website>`_ for making website changes.

New to Git/Github? Check out the `Newcomers' Guide <https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html>`_ and `Git Cheat Sheets <https://training.github.com/>`_.

Prefer a visual Git interface? Try `Github Desktop <https://desktop.github.com/>`_.

This guide is partially based on Astropy's `Development Workflow <https://docs.astropy.org/en/latest/development/workflow/development_workflow.html>`_.
