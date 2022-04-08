.. _dev_guide:
{% if not is_development %}
.. _newcomers:
.. _remote_data:
{% endif %}

*****************
Developer's Guide
*****************

{% if is_development %}

This article describes the guidelines to be followed by developers working on sunpy.
If you are thinking of contributing to sunpy please read the following carefully.

We currently recommend the :ref:`newcomers` as the place to start.
This goes over the basics and has links to useful tutorials on git.

.. toctree::
   :maxdepth: 2

    Getting started
    ^^^^^^^^^^^^^^^
    .. toctree::
      :maxdepth: 4

      contents/newcomers

    ---

    Standards
    ^^^^^^^^^
    .. toctree::
      :maxdepth: 1

      contents/code_standards
      contents/tests
      contents/documentation
      contents/example_gallery
      contents/pr_review_procedure
      contents/units_quantities

    ---

    Management
    ^^^^^^^^^^
    .. toctree::
      :maxdepth: 1

      contents/dependencies
      contents/funding
      contents/maintainer_workflow
      contents/ci_jobs

    ---

    API
    ^^^
    .. toctree::
      :maxdepth: 1

      contents/api
      contents/config
      contents/extending_fido
      contents/logger
      contents/new_objects
      contents/remote_data

{%else%}

Please go `here <https://docs.sunpy.org/en/latest/dev_guide/index.html>`__ for our up to date developer's guide.

{%endif%}
