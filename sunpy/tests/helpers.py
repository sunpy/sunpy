# -*- coding: utf-8 -*-
import warnings

import pytest

__all__ = ['skip_glymur', 'warnings_as_errors']

# SunPy's JPEG2000 capabilities rely on the glymur library.  First we check to
# make sure that glymur imports correctly before proceeding.
try:
    import glymur
except ImportError:
    SKIP_GLYMUR = True
else:
    SKIP_GLYMUR = False

# Skip ana tests if we are on Windows or we can't import the c extension.
skip_glymur = pytest.mark.skipif(SKIP_GLYMUR, reason="Glymur can not be imported")

@pytest.fixture
def warnings_as_errors(request):
    warnings.simplefilter('error')

    request.addfinalizer(lambda *args: warnings.resetwarnings())
