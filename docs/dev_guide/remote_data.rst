.. _remote_data:

************
Remote Files
************

Remote Data Manager
===================

There are functions which require files to be downloaded from the internet before they are executed.
If the files are changed on the remote server or if they are tampered with after the file is downloaded, the function may return erraneous results.

Remote Data Manager provides the developers a way to download the files easily and ensure that the files are intact when they are used.

Usage
-----

The manager has to be imported before it is used::

    from sunpy.data import manager



`~sunpy.data.data_manager.manager.DataManager.require` is a decorator which is used to inform the data manager that the function requires a file for it's execution.
`~sunpy.data.data_manager.manager.DataManager.get` function is used to access the file inside the function.
Suppose a function requires a file with url 'http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt' which has a SHA1 hash of '0f56b28dd53a99556254e66ba2c0401d567d0e94'.
The following example will show how this function can be implemented.::


    @manager.require('test_file',
                    ['http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt'],
                    '0f56b28dd53a99556254e66ba2c0401d567d0e94')
    def test_function():
        file_path = manager.get('test_file')
        # do something with file
        return file_path

Cache
=====

Remote files sometimes have to be cached.
This is to ensure that the files are not redownloaded frequently, thus saving users' disc space and well as internet bandwidth.
The expiry of the cache is configured by the user.

Usage
-----

The following example shows how cache can be used::

    from sunpy.data import cache

    def test_function():
        file_path = cache.download('http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt')
        return file_path


    test_function()  # The file will be downloaded with this call
    test_function()  # subsequent calls won't redownload the file


Testing
-------

A pytest fixture is provided for ease of mocking network requests when using cache.
The following example demonstates the usage of the fixture.::

    @pytest.fixture()
    def local_cache(sunpy_cache):
        sunpy_cache = sunpy_cache('sunpy.test_module.cache')
        sunpy_cache.add('http://example.com/test_file',
                        'test_data_path')

The above snippet creates a pytest fixture called `local_cache`. This fixture can be used in wherever the files have to be mocked.
An example is given below.::

    def test_test_function(local_cache):
        # inside this function the mocked cache is used

        # test_function uses 'http://example.com/test_file'
        assert test_function() == True
