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

        .. toctree::
            :maxdepth: 2

            getting_started

    .. grid-item-card::
        :class-card: card

        .. toctree::
            :maxdepth: 2

            conventions

    .. grid-item-card::
        :class-card: card

        .. toctree::
            :maxdepth: 2

            repo_management

    .. grid-item-card::
        :class-card: card

        .. toctree::
            :maxdepth: 2

            extending_sunpy

{%else%}

Please go `here <https://docs.sunpy.org/en/latest/dev_guide/index.html>`__ for our up to date developer's guide.

{%endif%}
