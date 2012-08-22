import weakref
import sunpy

class MapManager(weakref.WeakSet):
    """Weak referenced set of maps created using functions decorated by manage_maps."""
    pass

def manage_maps(fn):
    """Maps returned by functions decorated with manage_maps (eg. sunpy.make_map)
    will be registered in the sunpy.map_manager list."""

    def fn_manage_maps(*args, **kwargs):
        ret = fn(*args, **kwargs)
        sunpy.map_manager.add(ret)
        return ret
    return fn_manage_maps
