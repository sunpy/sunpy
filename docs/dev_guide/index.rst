.. _dev_guide:

{% if not is_development %}
.. _newcomers:
.. _remote_data:
{% endif %}

*****************
Developer's Guide
*****************

{% if is_development %}

This section contains the various guidelines to be followed by anyone working on sunpy.

.. grid:: 1 2 2 2
    :gutter: 3

    .. grid-item-card::
        :class-card: card

        Getting started
        ^^^^^^^^^^^^^^^

        * :doc:`contents/newcomers`
        * :doc:`contents/pr_checklist`
        * :doc:`contents/conda_for_dependencies`

    .. grid-item-card::
        :class-card: card

        Conventions
        ^^^^^^^^^^^

        * :doc:`contents/ai_usage`
        * :doc:`contents/code_standards`
        * :doc:`contents/tests`
        * :doc:`contents/documentation`
        * :doc:`contents/example_gallery`
        * :doc:`contents/pr_review_procedure`
        * :doc:`contents/units_quantities`

    .. grid-item-card::
        :class-card: card

        Repo management
        ^^^^^^^^^^^^^^^

        * :doc:`contents/maintainer_workflow`
        * :doc:`contents/dependencies`
        * :doc:`contents/ci_jobs`
        * :doc:`contents/backports`
        * :doc:`contents/deprecation`
        * :doc:`contents/funding`

    .. grid-item-card::
        :class-card: card

        Extending SunPy
        ^^^^^^^^^^^^^^^

        * :doc:`contents/public_api`
        * :doc:`contents/logger`
        * :doc:`contents/remote_data_manager`
        * :doc:`contents/extending_fido`
        * :doc:`contents/new_map_class`
        * :doc:`contents/custom_map_rotate`
        * :doc:`contents/scraper_migration`
        * :doc:`contents/soar/index`

.. toctree::
    :hidden:
    :maxdepth: 2

    getting_started
    conventions
    repo_management
    extending_sunpy

{%else%}

Please go `here <https://docs.sunpy.org/en/latest/dev_guide/index.html>`__ for our up to date developer's guide.

{%endif%}
