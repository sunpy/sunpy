.. _dev_guide_remote_data_manager_tests:

Tests using the Remote Data Manager
===================================

A pytest fixture (``sunpy_cache``) is provided for ease of mocking network requests when using cache.
The following example demonstrates the usage of the fixture:

.. code-block:: python

    @manager.require('test_file',
                     ['http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt'],
                     '4c85b04a5528aa97eb84a087450eda0421c71833820576330bba148564089b11')
    def test_function():
        return manager.get('test_file')

    @pytest.fixture()
    def local_cache(sunpy_cache):
        sunpy_cache = sunpy_cache('sunpy.test_module.cache')
        sunpy_cache.add('http://example.com/test_file',
                        'test_data_path')

The above snippet creates a pytest fixture called ``local_cache``.
This fixture can be used in wherever the files have to be mocked.
An example is given below:

.. code-block:: python

    def test_test_function(local_cache):
        # inside this function the mocked cache is used

        # test_function uses 'http://example.com/test_file'
        assert test_function() == True
