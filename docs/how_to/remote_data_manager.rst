.. _sunpy-how-to-the-remote-data-manager:

***************************
Use the Remote Data Manager
***************************

Often, data prep or analysis functions require files to be downloaded from a remote server.
The remote data manager handles the usage of remote data files including file verification using hashes.
For example, say a function, ``test_function``, requires the remote data file: `predicted-sunspot-radio-flux.txt <http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt>`__, which has the SHA256 hash "4c85b04a5528aa97eb84a087450eda0421c71833820576330bba148564089b11".
To ensure that this exact version of the file is downloaded when ``test_function`` is called, use the `~sunpy.data.data_manager.manager.DataManager.require` decorator.

.. code-block:: python

    >>> from sunpy.data import manager

    >>> @manager.require('test_file',
    ...                  ['http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt'],
    ...                  '4c85b04a5528aa97eb84a087450eda0421c71833820576330bba148564089b11')
    ... def test_function():
    ...     pass

The first time the function is called, the file will be downloaded and then cached such that subsequent function calls will not trigger a download.

.. code-block:: python

    >>> test_function()  # The file will be downloaded  # doctest: +REMOTE_DATA
    >>> test_function()  # No download  # doctest: +REMOTE_DATA

To access the downloaded file inside the function, use the :meth:`~sunpy.data.data_manager.manager.DataManager.get` method,

.. code-block:: python

    >>> @manager.require('test_file',
    ...                  ['http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt'],
    ...                  '4c85b04a5528aa97eb84a087450eda0421c71833820576330bba148564089b11')
    ... def test_function():
    ...     return manager.get('test_file')

To call a function that uses the data manager, but skip the hash check, use the `~sunpy.data.data_manager.manager.DataManager.skip_hash_check` context manager.

.. code-block:: python

    >>> with manager.skip_hash_check():
    ...     test_function()  # doctest: +REMOTE_DATA
    PosixPath('.../sunpy/data_manager/predicted-sunspot-radio-flux.txt')

To replace the required file with another version, use the `~sunpy.data.data_manager.manager.DataManager.override_file` context manager.

.. code-block:: python

    >>> with manager.override_file('test_file', 'http://data.sunpy.org/sample-data/AIA20110319_105400_0171.fits'):
    ...     test_function()  # doctest: +REMOTE_DATA
    PosixPath('.../sunpy/data_manager/AIA20110319_105400_0171.fits')
