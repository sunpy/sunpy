`sunpy.net.Scraper` now includes treats files as spanning a full interval equal to the smallest increment specified.
For example, a pattern like ``"%Y.txt"`` that only contains a year interval will be considered to span that full year.

This means searches that fall entirely within the interval will return that file where previously they did not.
For example matching ``"%Y.txt"`` with ``TimeRange('2022-02-01', '2022-04-01')`` will now return ``["2022.txt"]`` where previously no files were returned.
