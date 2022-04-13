.. _map_rotate_custom:

****************************
Adding a new rotation method
****************************

If you want to add a new rotation method that can be used by `sunpy.map.GenericMap.rotate`, you can do the following:

.. code-block:: python

    from sunpy.image.transform import add_rotation_function

    @add_rotation_function("my_rotate",
                        handles_clipping=False, handles_image_nans=False, handles_nan_missing=False)
    def _rotation_mine(image, matrix, shift, order, missing, clip):
        # Rotation code goes here
        return rotated_image

Then you can do:

.. code-block:: python

    from sunpy.map import Map
    from sunpy.data import sample

    hmi_map = sunpy.map.Map(sample.HMI_LOS_IMAGE)

    rot_map = hmi_map.rotate(order=3, recenter=True, method="my_rotate")

There are several things you will want to be aware of when adding a new rotation method.

* The `.add_rotation_function` means you can handle nans, clipping, and missing values in the image array
  before it is passed to your rotation function via the boolean flags.
  The documentation in `.add_rotation_function` explain how these work.
* The arguments into the rotation method must be ``image, matrix, shift, order, missing, clip``.
  There is no variation allowed.
* The return value of the rotation method must the rotated image array.
* You must describe the rotation method in the docstring of the function using unnumbered bullet points.

Below is a semi-functioning OpenCV2 rotation method as an example:

.. code-block:: python

    from sunpy.image.transform import add_rotation_function

    @add_rotation_function("cv2",
                        handles_clipping=False, handles_image_nans=False, handles_nan_missing=False)
    def _rotation_cv2(image, matrix, shift, order, missing, clip):
        """
        * Rotates using :func:`cv2.warpAffine`
        * This is only implemented for three orders of interpolation, 0, 1 and 3.
        * An input image with integer data is cast to floats prior to passing to
          :funccv2.warpAffine`. The output image can be re-cast using
          :meth:`numpy.ndarray.astype` if desired.
        * Does not let :func:`~cv2.warpAffine` handle clipping due to
          inconsistent handling across interpolation orders
        * Does not let :func:`~cv2.warpAffine` handle image NaNs because they
          are not supported.
        * Does not pass NaN as ``missing`` to :func:`~cv2.warpAffine` due to
          unsupported input.
        """
        try:
            import cv2
        except ImportError:
            raise ImportError("The opencv-python package is required to use this rotation method.")

        # OpenCV warpAffine does not support order 2,4,5
        _CV_ORDER_FLAGS = {
            0: cv2.INTER_NEAREST,
            1: cv2.INTER_LINEAR,
            3: cv2.INTER_CUBIC,
        }
        try:
            order = _CV_ORDER_FLAGS[order]
        except KeyError:
            raise ValueError(f"Input order={order} not supported in openCV. ",
                            "Please use order=0, 1, or 3.")
        trans = np.eye(3, 3)
        rot_scale = np.eye(3, 3)
        trans[:2, 2] = [-shift[0], -shift[1]]
        rot_scale[:2, :2] = matrix.T
        rmatrix = (rot_scale @ trans)[:2]
        if issubclass(image.dtype.type, numbers.Integral):
            warn_user("Integer input data has been cast to float64.")
            adjusted_image = image.astype(np.float64)
        else:
            adjusted_image = image.copy()
        h, w = adjusted_image.shape
        return cv2.warpAffine(adjusted_image, rmatrix, (w, h), flags=order,
                            borderMode=cv2.BORDER_CONSTANT, borderValue=missing)
