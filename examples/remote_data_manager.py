"""
========================================
Using sunpy.data.manager
========================================

This is an example to show how to use remote data manager to
handle remote data in sunpy.
"""
# Remote data manager is used to handle the usage of remote files in the code
# base in a sane way.
##############################################################################
from sunpy.data import manager

##############################################################################
# Let's start by defining a function which uses some remote files.
# Let's assume the function requires a file, say http://212.183.159.230/5MB.zip, which has
# the SHA1 hash '0cc897be1958c0f44371a8ff3dddbc092ff530d0'.

@manager.require('test_file',
                 ['http://212.183.159.230/5MB.zip'],
                 '0cc897be1958c0f44371a8ff3dddbc092ff530d0')
def test_function():
    pass


##############################################################################
# To access the downloaded file inside the function, you can use `manager.get` function
# `manager.get` returns a `pathlib.Path` object.

@manager.require('test_file',
                 ['http://212.183.159.230/5MB.zip'],
                 '0cc897be1958c0f44371a8ff3dddbc092ff530d0')
def test_function():
    return manager.get('test_file')


##############################################################################
# The first time the function is called, the file will be downloaded.
# During subsequent, no downloading will take place.

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

with manager.replace_file('test_file', 'http://212.183.159.230/10MB.zip'):
    test_function()
