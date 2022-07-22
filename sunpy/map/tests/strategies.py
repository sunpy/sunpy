import hypothesis.strategies as st
import numpy as np
from hypothesis import assume
from hypothesis.extra.numpy import arrays


@st.composite
def matrix_meta(draw, key):
    """
    Create an arbitrary but valid (ie non-singular) PCi_j or CDi_j matrix.

    Parameters
    ----------
    key : {'pc', 'cd'}
    """
    arr = draw(arrays(
        float, (2, 2),
        elements=st.floats(min_value=-1, max_value=1, allow_nan=False))
    )
    # Make sure matrix isn't singular by manually computing the determinant
    assume(np.abs(arr[1, 1]*arr[0, 0] - arr[0, 1]*arr[1, 0]) > 1e-8)
    return {f'{key}1_1': arr[0, 0],
            f'{key}1_2': arr[0, 1],
            f'{key}2_1': arr[1, 0],
            f'{key}2_2': arr[1, 1]}
