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

        .. toctree::
            :maxdepth: 3

            contents/newcomers

    .. grid-item-card::
        :class-card: card

        Conventions
        ^^^^^^^^^^^

        .. toctree::
            :maxdepth: 1

            contents/code_standards
            contents/tests
            contents/documentation
            contents/example_gallery
            contents/pr_review_procedure
            contents/units_quantities

    .. grid-item-card::
        :class-card: card

        Repo management
        ^^^^^^^^^^^^^^^

        .. toctree::
            :maxdepth: 1

            contents/maintainer_workflow
            contents/dependencies
            contents/ci_jobs
            contents/funding

    .. grid-item-card::
        :class-card: card

        Extending sunpy
        ^^^^^^^^^^^^^^^

        .. toctree::
            :maxdepth: 1

            contents/public_api
            contents/extending_fido
            contents/logger
            contents/new_objects
            contents/remote_data
            contents/map_rotate_custom

{%else%}

Please go `here <https://docs.sunpy.org/en/latest/dev_guide/index.html>`__ for our up to date developer's guide.

{%endif%}
