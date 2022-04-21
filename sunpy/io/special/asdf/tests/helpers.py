from io import BytesIO

import asdf

# Remove this file when we depend on asdf 2.9+


def roundtrip_object(obj, version=None):
    """
    Add the specified object to an AsdfFile's tree, write the file to
    a buffer, then read it back in and return the deserialized object.

    Parameters
    ----------
    obj : object
        Object to serialize.
    version : str or None.
        ASDF Standard version.  If None, use the library's default version.

    Returns
    -------
    object
        The deserialized object.
    """
    buff = BytesIO()
    with asdf.AsdfFile(version=version) as af:
        af["obj"] = obj
        af.write_to(buff)

    buff.seek(0)
    with asdf.open(buff, lazy_load=False, copy_arrays=True) as af:
        return af["obj"]
