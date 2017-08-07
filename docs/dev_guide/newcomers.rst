Newcomer's Guide
================

Welcome! If you came across this page, you just might be new to SunPy!

SunPy aims to be a comprehensive Python package that allows solar physicists to deep dive in the vast amount of solar data available.
The people who help develop or contribute to SunPy are varied in ability and experience with the vast majority being volunteers who dedicate time each week.
We pride ourselves on being a welcoming community and we would love to have you become a part of this!

If you have any questions, comments or just want to say hello, we have online chat on `Matrix`_ which requires no registration or a bit more old school, a `Google Group`_ which you can email.

.. _Matrix: https://riot.im/app/#/room/#sunpy-general:matrix.org
.. _Google Group: https://groups.google.com/forum/#!forum/sunpy

How To Contribute
-----------------

Not Code
########

Since SunPy is the work of many volunteers, we are always looking for more people to contribute however they can.
A common misconception (which applies to any package) is that all we really want is some Python code, in fact, we do not require only code contributions!
If you do not have the time or the desire to code, we have lots of areas where we can use help.

Most of the items within SunPy has been documented and is hosted `online here`_.
This documentation is not only for the code itself but contains setup instructions, quick start guides and worked examples.
However, documentation is never complete; there are always areas that could be expanded upon or could do with some proof reading to check whether what is written is easy to follow and understandable.
If parts are confusing or difficult to follow, we would love suggestions or improvements!

.. _online here: http://docs.sunpy.org/en/latest/dev.html

Code
####

If you would prefer to code Python instead, we have a SunPy issue list on Github `issues`_ where all the known issues with SunPy are kept.
An issue might have a label and the `Package Novice label`_ is a good place to start.
These are issues that have been deemed a good way to be eased into SunPy and are achievable with little understanding of the SunPy codebase.
If you came here from the Mozilla Sprint, we have a `mozsprint label.`_
These issues are similar in nature to the Package Novice label but with a slight difference that these issues are (we think) achievable within the timescale of a sprint, typically one day of hacking.

.. _issues: https://github.com/sunpy/sunpy/issues
.. _Package Novice label: https://github.com/sunpy/sunpy/issues?q=is%3Aissue+is%3Aopen+label%3Apackage-novice
.. _mozsprint label.: https://github.com/sunpy/sunpy/issues?q=is%3Aissue+is%3Aopen+label%3Amozsprint

How to setup a work environment
-------------------------------

We recommend using miniconda or Anaconda instead of your native operating system Python packages as it makes installing SunPy easier but also debugging.
If you have never setup a conda environment before, here is a very quick guide to getting setup.

The first step is to install the `miniconda`_ version that corresponds to your operating system.
Once you have downloaded the miniconda version you need, you can run the download.
If you get stuck, the install `instructions here`_ for each operating system.

Next we will want to setup the conda environment.
SunPy currently lives in the conda-forge channel which we can add this way.

.. code:: bash

    conda config --add channels conda-forge
    conda create -n sunpy-dev sunpy hypothesis pytest-mock pip
    source activate sunpy-dev

This will create a new conda environment called sunpy-dev and install the latest version of SunPy from the conda-forge channel.
The next step is remove the conda version of SunPy and install the in-development version of SunPy.
This requires that `git`_ be installed before hand.

Note that if you have a `GitHub`_ account, you can `fork`_ the `SunPy repository`_ (the fork button is to the top right) and use that url
for the ``git clone``.
This will make submitting changes easier in the long term for you.

.. code:: bash

    conda remove sunpy
    git clone https://github.com/sunpy/sunpy.git sunpy-git
    cd sunpy-git
    pip install -e .

Now you have the latest version of SunPy installed and are ready to work on it using your favorite editor!
Ideally, when you start making changes you want to create a git branch to work on:

.. code:: bash

    git checkout -b my_fix

You can change my\_fix to anything you prefer.
If you get stuck or want help, just `ask here`_!

.. _miniconda: https://conda.io/miniconda.html
.. _instructions here: https://conda.io/docs/install/quick.html#quick-install
.. _git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
.. _GitHub: https://github.com/
.. _fork: https://guides.github.com/activities/forking/
.. _SunPy repository: https://github.com/sunpy/sunpy
.. _ask here: https://riot.im/app/#/room/#sunpy-general:matrix.org

Send it back to us
------------------

Once you have some changes you would like to send to us.
To start you would need to commit the changes.

.. code:: bash

    git commit -a -m '<message>'

Where you replace ``<message>`` with some text of the work you have done.

Next step is to submit the changes back to SunPy.

The preferred method is that you submit a Pull Request (PR) using GitHub.
This will submit the code to SunPy where we can view the changes but also the inbuilt GitHub helpers allow some automatic review of the
submitted code.
If you are new to pull requests, here is a `friendly guide`_.
This way, we can review the code as a community and offer suggestions or accept it!

If you do not have time to finish what you started on or ran out of time during a sprint and do not want to submit a pull request, you can
create a git patch and send it to the `Google Group`_ or `email a SunPy contributor`_.
This way, you still get acknowledged for the work you did and this can be viewed within the SunPy git history.

.. code:: bash

    git format-patch master --stdout > my_fix.patch

You can rename ``my_fix`` to something more relevant to what you did.
This can be sent to a contributor or attached in the Google group.

Just remember, if you hit any problems get in touch!

Finally, a in-depth version of this guide is `located here`_.

.. _friendly guide: https://guides.github.com/activities/hello-world/
.. _Google Group: https://groups.google.com/forum/#!forum/sunpy
.. _email a SunPy contributor: stuart@mumford.me.uk
.. _located here: http://docs.sunpy.org/en/stable/dev.html
