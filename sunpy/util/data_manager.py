import weakref
import sunpy

class DataManager(weakref.WeakSet):
    """Weak referenced set of SunPy data types created using functions decorated by manage_data.
    The list is at sunpy._data_manager"""
    pass

def manage_data(fn):
    """SunPy data types returned by functions decorated with manage_data (eg. sunpy.make_map)
    will be registered in the sunpy._data_manager list."""

    def fn_manage_data(*args, **kwargs):
        ret = fn(*args, **kwargs)
        sunpy.data_manager.add(ret)
        return ret
    return fn_manage_data
