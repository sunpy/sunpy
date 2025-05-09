.. _newcomers:

****************
Newcomers' Guide
****************

If you have come across this page, you just might be new to the SunPy project - thanks for your interest in contributing to SunPy!
SunPy is an open project that encourages anyone to contribute in any way possible.
The people who help develop or contribute to SunPy are varied in ability and experience, with the vast majority being volunteers.

If you have any questions, comments, or just want to say hello, we have an active `chat room`_, a `mailing list`_ and a `community forum`_.

.. _chat room: https://app.element.io/#/room/#sunpy:openastronomy.org
.. _mailing list: https://groups.google.com/forum/#!forum/sunpy
.. _community forum: https://community.openastronomy.org/c/sunpy/5

How to Contribute to sunpy
==========================

.. _newcomers_generative_ai:

.. note:: Usage of Generative AI

    We expect authentic engagement in our community.
    Be wary of posting output from Large Language Models or similar generative AI as comments on GitHub or any other platform, as such comments tend to be formulaic and low quality content.
    If you use generative AI tools as an aid in developing code or documentation changes, ensure that you fully understand the proposed changes and can explain why they are the correct approach and an improvement to the current state.

Not Code
--------

A common misconception (which applies to any package) is that all we really want is some Python code, in fact, we do not require only code contributions!
If you do not have the time or the desire to code, we have severals of areas where we can use help.

Reporting Issues
^^^^^^^^^^^^^^^^

If you use sunpy and stumble upon a problem, the best way to report it is by opening an `issue`_ on our GitHub issue tracker.
This way we can help you work around the problem and hopefully fix the problem!

You will need to sign into `GitHub`_ to report an issue and if you are not already a member of Github, you will have to join.
Joining GitHub will make it easier to report and track issues in the future.

If you do not want to join Github, then another way to report your issue is email the SunPy developers list `sunpy-dev@googlegroups.com`_.

When reporting an issue, please try to provide a short description of the issue with a small code sample, this way we can attempt to reproduce the error.
Also provide any error output generated when you encountered the issue, we can use this information to debug the issue.
For a good example of how to do this see issue `#2879`_.

If there is functionality that is not currently available in sunpy you can make a feature request.
Please write as much information as possible regarding the feature you would like to see in sunpy.

When you go to open either an issue or a feature request, there is a GitHub template that will guide you on the type of information we are seeking.
Please be sure to read over the comments in the GitHub text box.

.. _issue: https://github.com/sunpy/sunpy/issues
.. _sunpy-dev@googlegroups.com: https://groups.google.com/forum/#!forum/sunpy-dev
.. _#2879: https://github.com/sunpy/sunpy/issues/2879

Documentation
^^^^^^^^^^^^^

sunpy has `online documentation`_ and we try to make sure its as comprehensive as possible.
This documentation contains the API of sunpy but also a user guide, an example gallery and developer documents.

However, documentation for any project is a living document.
It is never complete and there are always areas that could be expanded upon or could do with proof reading to check if the text is easy to follow and understandable.
If parts are confusing or difficult to follow, we would love suggestions or improvements!

.. _online documentation: https://docs.sunpy.org/en/latest/index.html

Reviewing a Pull Request
^^^^^^^^^^^^^^^^^^^^^^^^

We at any one time have a variety of `pull requests`_ open and getting reviews is important.
Generally the more people that can look over a pull request the better it will turn out and we encourage everyone to do so.

.. _pull requests: https://github.com/sunpy/sunpy/pulls

Code
----

If you would prefer to code instead, the best way to start is to work on an existing and known `issues`_.
We have several repositories you can investigate.
The main one is the sunpy repository with where all the known `issues`_ with sunpy are detailed.
Each issue should have a series of labels that provide information about the nature of the issue.
If you find an issue you'd like to work on, please make sure to add a comment to let people know that you are working on it! This will make it less likely that effort is duplicated.

.. note::

    sunpy is Free and open-source software (FOSS), under the BSD-2 license. By contributing you are stating that you have the right to and agree to have your work distributed under the terms of this license.

    This applies to everyone who wants to contribute during work time no matter who their employer is.
    You should start by checking if there is a Open Source Software Policy (e.g., `Stanford's policy <https://otl.stanford.edu/open-source-stanford>`__) for your work place.
    If not, `OSS-Watch <http://oss-watch.ac.uk/resources/contributing>`__ summarises what you will need to check and who to ask, however this resource is aimed at a UK readers.
    As an example, `Stanford's guidance <https://otl.stanford.edu/software>`__ allows someone to contribute and open source their code.
    If you are unsure if your university or institution allows you to contribute under the BSD-2 license, you should contact the relevant department or administrator that deals with copyright at your institution.

If you are unsure where to start we suggest the `Good First Issue label`_.
These are issues that have been deemed a good way to be eased into sunpy and are achievable with little understanding of the sunpy codebase.
Please be aware that this does not mean the issue is "easy", just that you do not need to be aware of the underlying structure of sunpy.

We also tag issues for specific events such as  `Hacktoberfest`_ under the `Hacktoberfest label`_.
The scope of the issues should be appropriate for that specific event.
We do participate in several other events but right now we do not have dedicated labels.
So please use the above labels for starting issues!

In addition, we have several other repositories that have open issues and you might find these more interesting than the main repository.

Python:

* `ndcube <https://github.com/sunpy/ndcube>`__
* `drms <https://github.com/sunpy/drms>`__
* `radiospectra <https://github.com/sunpy/radiospectra>`__
* `ablog <https://github.com/sunpy/ablog>`__
* `irispy <https://github.com/sunpy/irispy>`__
* `sunkit-image <https://github.com/sunpy/sunkit-image>`__

CSS/HTML:

* `sunpy-sphinx-theme <https://github.com/sunpy/sunpy-sphinx-theme>`__
* `sunpy.org <https://github.com/sunpy/sunpy.org>`__

.. _issues: https://github.com/sunpy/sunpy/issues
.. _Good First Issue label: https://github.com/sunpy/sunpy/issues?utf8=%E2%9C%93&q=is%3Aissue+is%3Aopen+label%3A%22Good+First+Issue%22
.. _Hacktoberfest: https://hacktoberfest.digitalocean.com/
.. _Hacktoberfest label: https://github.com/sunpy/sunpy/issues?q=is%3Aissue+is%3Aopen+label%3AHacktoberfest

If you already have code that you've developed or already have your own idea of code you'd like to work on please first have a look at the issue list to see if any existing issues are related to your idea.
If you find nothing then create your own issue to stimulate the addition of your code and make sure to let people know about it chat room or by email.
Creating an issue creates a permanent record.
It may be possible your idea may be outside of the scope of the repository you are trying to contribute to and the issue comments are a great place to have that discussion where potential future contributors can also see.

Setting up a development environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The instructions in the following section are based upon two resources:

* `Astropy Contributing Quickstart <https://docs.astropy.org/en/latest/development/quickstart.html>`__
* `Astropy Contribution Example <https://docs.astropy.org/en/latest/development/git_edit_workflow_examples.html>`__

**We strongly recommend that you read these links.**
These links are more in-depth than this guide but you will need to replace ``astropy`` with ``sunpy``.

In order to start coding you will need a local Python environment and we would recommend using `miniforge`_ (shortened to conda from here on).
This method will bypass your operating system Python packages and makes the entire process easier.

The first step is to install the version of miniforge which matches your system you can follow the instructions on :ref:`sunpy-tutorial-installing-miniforge` sunpy's install guide, but do not follow the "Installing sunpy" section.
Next we will want to setup the conda environment:

.. code:: bash

    # Note you might need to add python=<version> if a new release of Python has come out very recently.
    # Typically it will take around 3 months before we can support the latest version of Python.
    $ conda create -n sunpy-dev pip
    $ conda activate sunpy-dev

This will create a new conda environment called "sunpy-dev" and install the latest version of pip from the conda-forge channel.
The next step is get a development version of sunpy.
This will require that `git`_ be installed.
If you have a `GitHub`_ account, we suggest that you `fork`_ the `sunpy repository`_ (the fork button is to the top right) and **use that url for the clone step**.
This will make submitting changes easier in the long term for you:

.. warning::

    Do not clone the sunpy repository into ``$HOME/sunpy``. Depending on the operating system this location is used to store downloaded data files.
    This will cause conflicts later on, so the last argument (``sunpy-git``) on the ``git clone`` line will become the local folder name of the cloned repository.
    Otherwise you are free to clone to any other location.

.. code:: bash

    $ git clone https://github.com/<username>/sunpy.git sunpy-git
    $ cd sunpy-git
    # This adds the main sunpy repository as a remote called "upstream".
    # This will help you keep your fork up to date with the main sunpy repository.
    $ git remote add upstream https://github.com/sunpy/sunpy.git
    # This retrieves the tags from the main sunpy repository, which are used for determining the version of your fork.
    $ git fetch --tags upstream
    $ pip install -e ".[dev]"

.. note::
    If this does not work, it could be due to a missing C compiler (e.g., ``gcc`` or ``clang``) that is required to build sunpy at install.
    Getting the compiler either from your system package manager, XCode or conda-forge should address this.

Now you have the latest version of sunpy installed and are ready to work on it using your favorite editor!
Ideally, when you start making changes you want to create a git branch:

.. code:: bash

    $ git checkout -b my_fix

You can change ``my_fix`` to anything you prefer.
If you get stuck or want help, just `ask here`_!

.. _miniforge: https://github.com/conda-forge/miniforge?tab=readme-ov-file#download
.. _conda-forge: https://conda-forge.org/
.. _git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
.. _GitHub: https://github.com/
.. _fork: https://guides.github.com/activities/forking/
.. _sunpy repository: https://github.com/sunpy/sunpy
.. _ask here: https://app.element.io/#/room/#sunpy:openastronomy.org

Checking the code you have written
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now that you have written some code to address an issue.
You will need to check two things:

1. The changes you have made are correct, i.e., it fixes a bug or the feature works.
   This requires you to run the code either manually or by writing/running a test function.
   `pytest`_ is the framework we use for this.

2. The changes you have made follow the correct coding style.
   We follow the `PEP8`_ style for all Python code and depending on your setup, you can use a `linter program <https://realpython.com/python-code-quality/#how-to-improve-python-code-quality>`__ to check your code.
   For documentation, we follow the `numpydoc style <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`__.

We provide more more detail about our :ref:`test suite and how to write tests <testing>`, and how to :ref:`create and style documentation <docs_guidelines>`.

.. _pytest: https://docs.pytest.org/en/latest/

Send it back to us
^^^^^^^^^^^^^^^^^^
Once you have some changes you would like to submit, you will need to commit the changes.
This is a three stage process:

1. Use ``git status`` to see that the only changes locally are the right ones.
2. Use ``git add <path to file>`` to add the changes to ``git``.
3. Use ``git commit -m <message>`` to label those changes.
4. Use ``git push`` to update your fork (copy) of sunpy on GitHub.

Here you replace ``<message>`` with some text of the work you have done.
We strongly recommend having a good commit message and this `commit guide`_ is worth reading.

Next step is to open a pull request on GitHub.
If you are new to pull requests, here is the `GitHub guide`_ that is a detailed walkthrough.
Go to the "pull requests" tab on **your fork** and pressing the large green "New pull request" button.
Now on the right side from the box marked "compare" you can select your branch.
Do one final check to make sure the code changes look correct and then press the green "Create pull request" button.

When you open your pull request, we have a GitHub template that will guide you on what to write in the message box.
Please fill this in and title the pull request.
Now the final step is to press the green "Create pull request" button.

As soon as you do this, you will be greeted by a message from the "sunpy bot" as well as several continuous integration checks.
These are explained on our :ref:`Pull Request Review <pr_review>` page.
But what is important to know is that these run a series of tests to make sure that the changes do not cause any new errors.
We strongly recommend that any code changes you have had, follow the `PEP8`_ style and that you have ran the code locally to make sure any changes do not break any existing code.
We provide an overview on how to run the test suite :ref:`here <testing>`.
Now we (the sunpy community) can review the code and offer suggestions and once we are happy, we can merge in the pull request.

If you do not have time to finish what you started on or ran out of time during a sprint and do not want to submit a pull request, you can create a git patch instead:

.. code:: bash

    $ git format-patch main --stdout > my_fix.patch

You can rename ``my_fix`` to something more relevant.
This way, you still get acknowledged for the work you have achieved.
Now you can email this patch to the  `Google Group`_ .

Just remember, if you have any problems get in touch!

.. _commit guide: https://chris.beams.io/posts/git-commit/
.. _GitHub guide: https://guides.github.com/activities/hello-world/
.. _PEP8: https://realpython.com/python-pep8/
.. _Google Group: https://groups.google.com/forum/#!forum/sunpy

Summer of Code(s)
^^^^^^^^^^^^^^^^^

If you are interested in a "Summer of Code" project with sunpy, we have information on our `wiki`_ which has guidelines, advice, application templates and more!
Our projects are located on our umbrella's organization website, `OpenAstronomy`_.

.. _wiki: https://github.com/sunpy/sunpy/wiki#summer-of-codes
.. _OpenAstronomy: https://openastronomy.org/gsoc/
