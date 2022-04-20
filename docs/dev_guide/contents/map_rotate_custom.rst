.. _map_rotate_custom:

****************************
Adding a new rotation method
****************************

If you want to add a new rotation method that can be used by `sunpy.map.GenericMap.rotate`, you can do the following:

.. code-block:: python

    from sunpy.image.transform import add_rotation_function

    @add_rotation_function("my_rotate", handles_clipping=False, handles_image_nans=False, handles_nan_missing=False)
    def _rotation_mine(image, matrix, shift, order, missing, clip):
        # Rotation code goes here
        return rotated_image

The best explanation of the arguments is found from the detailed function documentation:

.. autofunction:: sunpy.image.transform.add_rotation_function
   :noindex:

Then you can do:

.. code-block:: python

    from sunpy.map import Map
    from sunpy.data import sample

    hmi_map = sunpy.map.Map(sample.HMI_LOS_IMAGE)

    rot_map = hmi_map.rotate(order=3, recenter=True, method="my_rotate")

It is important to know what keywords from `sunpy.map.GenericMap.rotate` are used by the underlying rotation function.

* ``angle`` - This is passed to `sunpy.image.transform.affine_transform` as a rotation matrix.
* ``rmatrix`` - This is passed to `sunpy.image.transform.affine_transform`
* ``order`` - This is passed to `sunpy.image.transform.affine_transform`
* ``scale`` - This is passed to `sunpy.image.transform.affine_transform`
* ``recenter`` - This is passed to `sunpy.image.transform.affine_transform`
* ``missing`` - This is passed to `sunpy.image.transform.affine_transform`
* ``use_scipy`` - This is passed to `sunpy.image.transform.affine_transform`, it is deprecated and removal is planned for sunpy 4.1, you should use ``method`` instead.
* ``method`` - This is passed to `sunpy.image.transform.affine_transform`
* ``clip`` - This is passed to `sunpy.image.transform.affine_transform`

The methods that ``sunpy`` use are defined in `sunpy/image/transform.py <https://github.com/sunpy/sunpy/blob/main/sunpy/image/transform.py>`__`.
