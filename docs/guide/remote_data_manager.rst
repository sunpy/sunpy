.. _how_to_use_the_remote_data_manager:

Use the Remote Data Manager
===========================

Often, data prep or analysis functions require files to be downloaded from a remote server.
The remote data manager handles the usage of remote data files including file verification using hashes.
For example, say a function, ``test_function``, requires the remote data file http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt, which has the SHA256 hash '4c85b04a5528aa97eb84a087450eda0421c71833820576330bba148564089b11'.
To ensure that this exact version of the file is downloaded when ``test_function`` is called, use the `~sunpy.data.data_manager.manager.DataManager.require` decorator.

The remote data manager provides developers with a way to download the files easily and ensure that the files are intact when they are used.
If the file is changed on the remote server, the data manager will raise an error so that the user is made aware that the intended file is not available.
If the file has changed on disk, the data manager will redownload the file after warning the user.

    >>> from sunpy.data import manager
    >>> @manager.require('test_file',
    ...                  ['http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt'],
    ...                  '4c85b04a5528aa97eb84a087450eda0421c71833820576330bba148564089b11')
    ... def test_function():
    ...     pass

The first time the function is called, the file will be downloaded and then cached such that subsequent function calls will not trigger a download.

The manager has to be imported before it is used:

.. code-block:: python

    >>> test_function()  # The file will be downloaded  # doctest: +REMOTE_DATA
    >>> test_function()  # No downloading here  # doctest: +REMOTE_DATA

`~sunpy.data.data_manager.manager.DataManager.require` is a decorator which is used to inform the data manager that the function requires a file for its execution.
`~sunpy.data.data_manager.manager.DataManager.get` function is used to access the file inside the function.
Suppose a function requires a file with url 'http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt' which has a SHA256 hash of '4c85b04a5528aa97eb84a087450eda0421c71833820576330bba148564089b11'.
The following example will show how this function can be implemented:

.. code-block:: python

.. code-block:: python

    >>> @manager.require('test_file',
    ...                  ['http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt'],
    ...                  '4c85b04a5528aa97eb84a087450eda0421c71833820576330bba148564089b11')
    ... def test_function():
    ...     return manager.get('test_file')

To call a function that uses the data manager, but skip the hash check, use the `~sunpy.data.data_manager.manager.DataManager.skip_hash_check` context manager.

.. code-block:: python

The following example shows how cache can be used:

.. code-block:: python

To replace the required file with another version, use the `~sunpy.data.data_manager.manager.DataManager.override_file` context manager.

.. code-block:: python


    test_function()  # The file will be downloaded with this call
    test_function()  # subsequent calls won't redownload the file

Testing
=======

A pytest fixture is provided for ease of mocking network requests when using cache.
The following example demonstrates the usage of the fixture:

.. code-block:: python

    @pytest.fixture()
    def local_cache(sunpy_cache):
        sunpy_cache = sunpy_cache('sunpy.test_module.cache')
        sunpy_cache.add('http://example.com/test_file',
                        'test_data_path')

The above snippet creates a pytest fixture called ``local_cache``. This fixture can be used in wherever the files have to be mocked.
An example is given below:

.. code-block:: python

    def test_test_function(local_cache):
        # inside this function the mocked cache is used

        # test_function uses 'http://example.com/test_file'
        assert test_function() == True
