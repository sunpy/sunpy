.. _version_control:

Version Control
===============

Source-code for SunPy is managed using `Git <https://git-scm.com>`_,
a Distributed Version Control system. Code branches are hosted on
`GitHub.com <https://github.com/sunpy/sunpy>`_, a free project hosting website
for Open-Source software. To contribute to SunPy you will need to create an
account on GitHub.

Creating Your Own Repo
----------------------

Each person contributing to SunPy should create their own code repository (repo)
on GitHub by forking the master repository or repo. All development is then
done on your fork, using topic branches to isolate work on different features.
New contributors can then initiate pull requests to have their code incorporated
into the SunPy master repository. Regular contributors may become members of
the `SunPy team <https://github.com/sunpy>`_ on GitHub. Code will be reviewed
by and comments will be provided to improve the code before it is accepted.

Getting Started
###############

Creating your own fork of SunPy on GitHub is easy to do. Go ahead and create
a free account on `create an account on GitHub <https://github.com/join>`_.
Github has some `great resources to help <https://help.github.com/>`_.
Here is a quick overview of the process.

Adding an SSH key to GitHub
###########################

Next, you need to tell GitHub who you are. In order to push any code to GitHub
you need to create a public SSH key and associate it with your GitHub account.
For instructions on how this is done, see the article on GitHub on
`Setting up git <https://help.github.com/set-up-git-redirect>`_ under
"Set Up SSH Keys". You only need to do this once, although if you plan to
work from multiple computers you will need to go through the process for each
computer.

Using HTTPS
###########

If you do not fancy using SSH you can access GitHub using HTTP/HTTPS.
A few things to note. Using HTTP only allows cloning of public repositories,
while HTTPS allows cloning of private repositories but also allows you to have
push access. This way you can type in your username and password to access
your repositories.

Identifying yourself
####################

Begin by identifying yourself to git (so all of your commits are signed with
this information): ::

 git config --global user.name "Firstname Lastname"
 git config --global user.email "your_email@youremail.com"

Forking SunPy
#############

Each contributor to SunPy must have their own copy of the SunPy master repo.
When working on the code, changes are made to your copy, and only when the
changes are completed, and have been verified to work, should they be
requested to be merged into the SunPy code base (through a pull request).
GitHub provides a simple mechanism to setup your own
personal repo by providing an option to `fork a repository
<https://help.github.com/fork-a-repo/>`_. When you create a fork of a GitHub
project, a copy of the repo will automatically be created for you, and a link
will be provided which you can use to download the code to your machine and
begin working on it.

To begin, fork the main SunPy repo on GitHub by clicking on the `Fork` button
on the `SunPy project page <https://github.com/sunpy/sunpy>`_

Next, you need to download the forked repository. Clone the fork to your
local machine, edit and run: ::

 git clone git@github.com:<your_username>/sunpy.git

or: ::

 git clone https://github.com/sunpy/sunpy.git

By default your fork of the repo on GitHub is identified by the name `origin`.
In order to keep the fork up to date with the main repo, it is useful to add it
as a `remote` in git: ::

 git remote add upstream https://github.com/sunpy/sunpy.git

To stay up to date regularly grab the latest changes to the SunPy master using
the commands: ::

 git fetch upstream master
 git merge upstream/master

This will merge the upstream code with your code so you don't
need to worry about it overwriting your changes. Next let's test this by making
a change to SunPy. First (always do this!) create a new branch to contain your
change.

 git checkout -b test

Go ahead and modify one of the files, or create a new file
(and then run :command:`git add`).

Commit and push the changes to your GitHub: ::

 git commit -a -m "My first commit"
 git push

You local repo is now synced with GitHub and ahead of the main repo as it
contains your change. Remember to commit after you've done a unit of work (i.e.
often). This will make it easier for you (in the future) and everyone else to
understand what you are doing. Also make sure to make your commit statements
clear and understandable to yourself and others. You can now go back to your
unchanged version of SunPy on your master branch by running the command

 git checkout master

Remember to always create a new branch to contain a certain set of work.

Installing SunPy
################

In order to use the version of SunPy located in your personal repository.
You need to install it using the `setup.py` script located in the top-level
folder. The `setup.py` script has several flags: ::
`develop` : Installs SunPy and builds all external libraries.
`build` or `build_ext`:  (Re)Builds the external libraries.
`clean --all`: Cleans all build files

Use the `setup.py` script like so: ::

 pip install -e ./

If you are interested in having different versions of sunpy in your
machine and you want to switch from one to another you could use
virtual environments. This is an easy task if you used conda as your
python package manager.

After a standard conda installation, :ref:`assuming you have also installed
the latest stable version of sunpy <main-install>`, you then proceed to create a new environment
as::

 conda create -n sunpy-dev python=3 sunpy

This will create a new environment called `sunpy-dev` with all of the
dependencies needed by sunpy. We then proceed to change to the new
environment::

 source activate sunpy-dev

Then we need to remove the stable version from this environment ::

 conda remove sunpy

to then install the version in your git repository ::

 cd to/sunpy/git/repository
 pip install -e ./

At this stage you can use the development version in which you are
working on. If you want to go back to the stable installation you can just
change the environment by ::

 source deactivate

Conclusion
##########

That's it! You now have your own personal SunPy repo to develop on. You could
hack away at it to your heart's content, pushing changes to your fork on
GitHub to share with others and to ensure that you have a backup online.

But what about when you want to start contributing back to the main SunPy
repo? That is the topic of the next section.

Branches
--------

Developers should create topic branches within their repos for most of their
main coding. Every repo starts with a single branch called `master`, which
seldom needs to be used. Instead, work on any particular feature, bug, or
portion of the code is done in its own separate branch. This way changes on
any particular issue are isolated from other unrelated changes. Users can even
work on several different branches simultaneously.

To create a new branch run: ::

 git branch branchname

To switch to the new branch: ::

 git checkout branchname

(or alternatively, :command:`git checkout -b branchname` will accomplish
the above).

Developers should create new branches for the features they are working on.
When they have finished making changes and the code has been tested and
verified to be working well, the code can be merged back into the SunPy
repo. This is usually done through something called a pull request.

Example Workflow
################

Here is an example workflow for a SunPy developer on any given day. Before
beginning this tutorial, follow the above instructions to grab a copy of the
SunPy repo.

Grabbing other people's changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first thing you want to do before you start coding anything new is to pull
in the latest code that others have written since you last did any coding. To
do this, run `git pull`: ::

    git pull upstream master

This will ensure that you don't edit a file that has changed since your last pull
which will lead to merge conflicts later on.

Code away
^^^^^^^^^

Assuming there are no merge conflicts (which shouldn't happen unless two people
are working on the same part of the same file), then you are ready to begin
coding. If there are conflicts check out our conflicts section.

Push your changes to GitHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^

As you code away on your local repo, you will need to keep git aware of what you are doing
and also your remote copy up to date.

To add a file, create the file then run: ::

    git add <yourfilename>

If you delete a file run: ::

    git rm <yourfilename>

To move a file: ::

    git mv <source> <destination>

To check to see if git is happy run: ::

    git status

which will give you a report of what has happened so far. Once you are at a good stopping point you should
"commit" your changes. This will provide you an opportunity to describe what you have done so far. To do this type: ::

    git commit -a -m "description of your changes"

After doing this you are ready to push your changes to your repo online with the command: ::

    git push

The local and remote copies of your repo are now synced.

Contributing to the main repo
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have made your desired changes, and committed and pushed your personal
branch, you need to decide whether or not to merge those changes back into the
main SunPy repo. If the changes you made are finished and have been tested and proven
stable (see the testing section below), then they can be merged into SunPy.
For now, lets assume that
your changes are complete and they are ready to be added to the main SunPy repo.
All contributed code to SunPy must be submitted as a "pull request". To do this go to the github
website and to your repo (remember to select the branch) then click on the "Pull
Request" button (in the upper right hand corner next to the Fork button which you've
used before). All initial pull requests must be made to the master branch unless they are a fix for specific version.
This will submit your code to a review. You will likely
receive some constructive comments on your code. To address these you can simply work
on your code and push those changes to your local repo. Those changes will be reflected
in your pull request. Once a member of
the SunPy dev team approves your pull request then your code will be
merged into the main SunPy repo
and your code will be part of the main SunPy code. Congratulations!

And that's it! It may seem like a lot at first but once you go through the
motions a few times it becomes very quick.

Conflict resolution
^^^^^^^^^^^^^^^^^^^

It may so happen that when you try to sync with the main repo there is a conflict error.
This means that someone else has been working on the same section of code
that you have. In such cases, the merge
command will issue a conflict warning and will then expect you do the merge
yourself. You can type: ::

   git mergetool

to go through the conflicts. This command will likely open some merging tools
which are already available on your computer. For example, on Mac OS X, it will open
FileMerge (if you have XCode installed). You can check on your progress by typing: ::

   git status

Once you are done, you should then commit your changes, in this case
the resolution of the conflict with: ::

   git commit -m "Resolved conflict between my and online version of file.py"

You can then proceed to push this change up to your branch.

Rebasing
^^^^^^^^

Sometimes it might be better to instead of merging in upstream/master, to rebase on top of upstream/master, or if you would like to clean up  your commit history if you deem it messy.
**However**, be warned that rebasing is a nuclear option.
If it goes wrong, it fundamentally changes your git history, there is no way back if you have not got a copy somewhere else, say your online fork of SunPy.
You can back out of a rebase during the process.

We will have a brief example here but since rebasing is a major step (depending on the complexity of the pull request) we would recommend checking out `this <https://www.digitalocean.com/community/tutorials/how-to-rebase-and-update-a-pull-request>`_ or `this one <https://www.atlassian.com/git/tutorials/rewriting-history/git-rebase>`_ tutorial.

If you are on your own branch and you have upstream added as a remote.
You can do ::

    git rebase upstream/master

which will rebase your commits on top of upstream/master and if there are no major changes, it should complete with no problem.
If you add a `-i`, this will turn on interactive mode ::

    git rebase -i upstream/master

you should see something like this ::

    pick 2231360 some old commit
    pick g3s62dc some mid commit you want to remove
    pick ee2adc2 Adds new feature
    # Rebase 2cf755d..ee2adc2 onto 2cf755d (9 commands)
    #
    # Commands:
    # p, pick = use commit
    # r, reword = use commit, but edit the commit message
    # e, edit = use commit, but stop for amending
    # s, squash = use commit, but meld into previous commit
    # f, fixup = like "squash", but discard this commit's log message
    # x, exec = run command (the rest of the line) using shell
    # d, drop = remove commit

Here you can change `pick` to any of the other commands that are listed and have that change the commits in your local history.
So if you wanted to remove the middle commit you would change ::

    pick g3s62dc some mid commit you want to remove

to ::

    drop g3s62dc some mid commit you want to remove

or if you wanted to keep the changes merge that commit into the previous commit ::

    squash g3s62dc some mid commit you want to remove

Now when you exit the screen, git will now apply the changes you are after.

If any problem arises, git will tell you and allow you either work through the problem using `git mergetool` or to abort the process `git rebase --abort`.

Backporting contribution
^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes a contribution needs to be backported to the latest stable branch, this
may be due to a bug being fixed or something similar.
There are different ways to do so, if the contribution contains just a couple
of commits, then the easiest is to `cherry-pick` them.
Assuming you are in the branch of your new feature (eg. `new_feature`), this
is what you need to do:

First you need to find out which commits you want to copy to the other branch: ::

  git log

Download/update the upstream branches to your local machine: ::

  git fetch upstream

Create a new branch from the version you want to backport, X.y: ::

  git checkout -b new_feature_X.y upstream/X.y

Copy the commits using `cherry-pick`, `xxxxxxxx` (`yyyyyyyy`) refers to the
oldest (newest) commit you want to backport. `^` at the end of the oldest is
to include it, otherwise will take the ones after that point: ::

  git cherry-pick xxxxxxxx^..yyyyyyyy

Push that new branch to your repository on github: ::

  git push origin new_feature_X.y

Once done, then you can create a new pull request to the X.y branch.
Remember to keep the same title that the original but adding [X.y] at the beginning.
Also add a reference to the original pull request in the comments with
the appropriate format: `#pr-number`.
