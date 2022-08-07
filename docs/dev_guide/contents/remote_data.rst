.. _remote_data:

*******************
Remote Data Manager
*******************

There are functions which require files to be downloaded from the internet before they are executed.
If the files change on the remote server or if they are tampered with after the file is downloaded, the function may return erroneous results.

The remote data manager provides developers with a way to download the files easily and ensure that the files are intact when they are used.
If the file is changed on the remote server, the data manager will raise an error so that the user is made aware that the intended file is not available.
If the file has changed on disk, data manager will re-download the file after warning the user.

Also there is a `~sunpy.util.hash_file` function.

Usage
=====

The manager has to be imported before it is used:

.. code-block:: python

    from sunpy.data import manager

`~sunpy.data.data_manager.manager.DataManager.require` is a decorator which is used to inform the data manager that the function requires a file for its execution.
`~sunpy.data.data_manager.manager.DataManager.get` function is used to access the file inside the function.
Suppose a function requires a file with url 'http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt' which has a SHA256 hash of '4c85b04a5528aa97eb84a087450eda0421c71833820576330bba148564089b11'.
The following example will show how this function can be implemented:

.. code-block:: python

    @manager.require('test_file',
                    ['http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt'],
                    '4c85b04a5528aa97eb84a087450eda0421c71833820576330bba148564089b11')
    def test_function():
        file_path = manager.get('test_file')
        # do something with file
        return file_path

Cache
-----

Remote files sometimes have to be cached.
This is to ensure that the files are not re-downloaded frequently, thus saving users' disk space and well as internet bandwidth.
The cache has an expiry time which is set by the ``sunpyrc`` config file.
To change this please see :ref:`customizing-sunpy`.

Usage
=====

The following example shows how cache can be used:

.. code-block:: python

    from sunpy.data import cache

    def test_function():
        file_path = cache.download('http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt')
        return file_path

    test_function()  # The file will be downloaded with this call
    test_function()  # subsequent calls won't re-download the file


Testing
=======

A **pytest** fixture is provided for ease of mocking network requests when using cache.
The following example demonstrates the usage of the fixture:

.. code-block:: python

    @pytest.fixture
    def local_cache(sunpy_cache):
        sunpy_cache = sunpy_cache('sunpy.test_module.cache')
        sunpy_cache.add('http://example.com/test_file', 'test_data_path')

The above snippet creates a pytest fixture called ``local_cache``.
This fixture can be used in wherever the files have to be mocked.
An example is given below:

.. code-block:: python

    def test_test_function(local_cache):
        # inside this function the mocked cache is used

        # test_function uses 'http://example.com/test_file'
        assert test_function()
