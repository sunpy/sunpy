.. _sunpy-topic-guide-add-a-new-rotation-method-to-map:

********************************
Add a new rotation method to Map
********************************

It is possible to select from a number of rotation methods when using :func:`sunpy.image.transform.affine_transform` and :meth:`sunpy.map.GenericMap.rotate`.
You can add a custom rotation method using the decorator :func:`sunpy.image.transform.add_rotation_function`:

.. code-block:: python

    from sunpy.image.transform import add_rotation_function

    @add_rotation_function("my_rotate", allowed_orders={0, 1, 3},
                           handles_clipping=False, handles_image_nans=False, handles_nan_missing=False)
    def _rotation_mine(image, matrix, shift, order, missing, clip):
        # Rotation code goes here
        return rotated_image

See the docstring for :func:`~sunpy.image.transform.add_rotation_function` for a detailed explanation of each of the decorator parameters and each of the required input parameters to the rotation function

Then you can do:

.. code-block:: python

    >>> from sunpy.map import Map
    >>> from sunpy.data import sample

    >>> hmi_map = sunpy.map.Map(sample.HMI_LOS_IMAGE)  # doctest: +SKIP
    >>> rot_map = hmi_map.rotate(order=3, recenter=True, method="my_rotate")  # doctest: +SKIP

The available rotation methods are all implemented using the :func:`~sunpy.image.transform.add_rotation_function` decorator, so you can look in `sunpy/image/transform.py <https://github.com/sunpy/sunpy/blob/main/sunpy/image/transform.py>`__ for examples of how to use this decorator.
