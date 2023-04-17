************************
sunpy core Documentation
************************

``sunpy`` is a community-developed, free and open-source solar data analysis environment for Python.
It provides search and download functionality, data containers for image and time series data, as well as commonly used coordinate frames and transformations between such frames.

.. image:: tour.svg
  :align: center

New Users
==========

For help installing ``sunpy`` refer to our :ref:`installation guide <installing>`.

.. todo::

    if #6615 is completed link here
    Introduction to sunpy core

If you are a first time user or new to Python, our **tutorial** provides a walkthrough on how to use the key features of sunpy core.
Once you've installed sunpy, start here.

* :ref:`tutorial`
    #. :ref:`installing`
    #. :ref:`units-sunpy`
    #. :ref:`time-in-sunpy`
    #. :ref:`coordinates-sunpy`
    #. :ref:`acquiring_data`
    #. :ref:`map_guide`
    #. :ref:`timeseries_guide`

Next Steps
==========

* :doc:`Example gallery <generated/gallery/index>` contains short-form guides on accomplishing common tasks using sunpy.
* :ref:`how_to_guide` contains short snippets of code for accomplishing specific tasks with sunpy. This section is most useful for answering "How do I..." questions.
* :ref:`guide` contains in-depth explanations of key topics and concepts in sunpy. This section of the documentation is more advanced than the others and is most useful for answering "why" questions.
* :ref:`API reference<reference>` contains a technical description of the inputs, outputs, and behavior of each component of sunpy core.

Further Info
============

* Find out :ref:`what's new <whatsnew>` with sunpy.
* View a list of :ref:`known issues <known_issues>` or report new issues on our `GitHub issue tracker`_.
* Join our `chat room`_ or post to our `discourse`_ for help.
* Read the :ref:`Developers guide <dev_guide>` to start contributing to the SunPy project.


.. _chat room: https://app.element.io/#/room/#sunpy:openastronomy.org
.. _GitHub issue tracker: https://github.com/sunpy/sunpy/issues
.. _discourse: https://community.openastronomy.org/c/sunpy/5

.. toctree::
    :maxdepth: 1
    :hidden:

    tutorial/index
    tutorial/installation
    generated/gallery/index
    how_to/index
    guide/index
    reference/index
    whatsnew/index
    about
    reference/known_issues
    reference/stability
    dev_guide/index
