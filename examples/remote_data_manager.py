"""
=========================
Using Remote Data Manager
=========================

This is an example to show how to use remote data manager to
handle remote data in Sunpy.
"""
##############################################################################
# Remote data manager is used to handle the usage of remote files in the code
#  base with file verification using hashes.

from sunpy.data import manager

##############################################################################
# Let's start by defining a function which uses some remote files.
# Let's assume the function requires a file,
# say http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt, which has
# the SHA1 hash '0f56b28dd53a99556254e66ba2c0401d567d0e94'.


@manager.require('test_file',
                 ['http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt'],
                 '0f56b28dd53a99556254e66ba2c0401d567d0e94')
def test_function():
    pass

##############################################################################
# To access the downloaded file inside the function, you can use
# `~sunpy.data.data_manager.manager.DataManager.get` function
# `manager.get` returns a `pathlib.Path` object.


@manager.require('test_file',
                 ['http://data.sunpy.org/sample-data/predicted-sunspot-radio-flux.txt'],
                 '0f56b28dd53a99556254e66ba2c0401d567d0e94')
def test_function():
    return manager.get('test_file')

##############################################################################
# The first time the function is called, the file will be downloaded.
# During subsequent calls, no downloading will take place.


test_function()  # The file will be downloaded
test_function()  # No downloading here

##############################################################################
# In case the user wants to skip the hash check, there is a helper context manager
# `~sunpy.data.data_manager.manager.DataManager.skip_hash_check`.

with manager.skip_hash_check():
    test_function()

##############################################################################
# If the user knows the function is going to use a file and want to replace it with another version
# they can do that too.

with manager.override_file('test_file', 'http://data.sunpy.org/sample-data/AIA20110319_105400_0171.fits'):
    test_function()
