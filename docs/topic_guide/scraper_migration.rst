.. _sunpy-topic-guide-scraper-migration:

************************************************************
Migrating to the new Scraper pattern system in ``sunpy`` 6.1
************************************************************

The `~sunpy.net.scraper.Scraper` has undergone a major change in the sunpy 6.1 release.

The old pattern system was regex-based, while the new pattern system is based on the `parse <https://github.com/r1chardj0n3s/parse/>`__
We hope this new system will be more flexible and easier to use.
`The readme of parse library provides a good overview of the syntax. <https://github.com/r1chardj0n3s/parse/#format-syntax>`__

How to migrate
==============

Previously, when defining a new Fido client or working with the `~sunpy.net.scraper.Scraper`, a regex-based ``baseurl`` and a parse-styled ``pattern`` attributes were required.
The new pattern system requires a single pattern string that covers the entire URL, a combination of the ``baseurl`` and ``pattern`` attributes from before.
This pattern string has to be passed as a ``format`` keyword to the `~sunpy.net.scraper.Scraper`.
During the deprecation period, the old ``pattern`` argument will still be accepted, but it is recommended to switch to the new pattern system.

Writing the new pattern
-----------------------

1. The new pattern covers the entire URL from the protocol, domain to the filename.
2. Instead of conveying time or numerical information in the URL in a datetime or regex format, it is now done in the `parse <https://github.com/r1chardj0n3s/parse/#format-syntax>`__ format.
   E.g., a ``"%Y"`` is now to be replaced with ``{{year:4d}}``.
   The keywords corresponding to the datetime format (the supported time keys - to be used within ``{{}}``) are:

   * 'year:4d'
   * 'year:2d'
   * 'month:2d'
   * 'month_name:l'
   * 'month_name_abbr:l'
   * 'day:2d'
   * 'day_of_year:3d'
   * 'hour:2d'
   * 'minute:2d'
   * 'second:2d'
   * 'microsecond:6d'
   * 'millisecond:3d'
   * 'week_number:2d'

3. The metadata attributes for extraction are written within double curly-braces ``{{}}`` for eg. ``{{ADAPTRealizations:3d}}``.
4. Single curly-braces ``{}`` are used in case of regular placeholders for Python format strings, with their respective values being passed as ``kwargs`` to `~sunpy.net.scraper.Scraper`.

Example
-------

Before migration:

.. code-block:: python

    baseurl = r'https://gong.nso.edu/adapt/maps/gong/%Y/adapt(\d){5}_(\d){2}(\w){1}(\d){3}_(\d){12}_(\w){1}(\d){8}(\w){1}(\d){1}\.fts\.gz'
    pattern = '{}adapt{ADAPTFileType:1d}{ADAPTLonType:1d}{ADAPTInputSource:1d}{ADAPTDataAssimilation:1d}{ADAPTResolution:1d}' + \
    '_{ADAPTVersionYear:2d}{ADAPTVersionMonth:1l}{ADAPTRealizations:3d}_{year:4d}{month:2d}{day:2d}{hour:2d}{minute:2d}' + \
    '_{ADAPTEvolutionMode:1l}{days_since_last_obs:2d}{hours_since_last_obs:2d}{minutes_since_last_obs:2d}{seconds_since_last_obs:2d}{ADAPTHelioData:1l}{ADAPTMagData:1d}.fts.gz'

After migration:

.. code-block:: python

    pattern = r'https://gong.nso.edu/adapt/maps/gong/{{year:4d}}/adapt{{ADAPTFileType:1d}}{{ADAPTLonType:1d}}{{ADAPTInputSource:1d}}{{ADAPTDataAssimilation:1d}}{{ADAPTResolution:1d}}' + \
    '_{{ADAPTVersionYear:2d}}{{ADAPTVersionMonth:1l}}{{ADAPTRealizations:3d}}_{{year:4d}}{{month:2d}}{{day:2d}}{{hour:2d}}{{minute:2d}}' + \
    '_{{ADAPTEvolutionMode:1l}}{{days_since_last_obs:2d}}{{hours_since_last_obs:2d}}{{minutes_since_last_obs:2d}}{{seconds_since_last_obs:2d}}{{ADAPTHelioData:1l}}{{ADAPTMagData:1d}}.fts.gz'
