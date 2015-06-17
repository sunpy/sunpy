from __future__ import unicode_literals
def get_package_data():
    return {'sunpy.data': ['sunpyrc'],
            'sunpy.data.test': ['*.*', '*/*.*']}
