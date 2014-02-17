# -*- coding: utf-8 -*-
import warnings

import pytest

@pytest.fixture
def warnings_as_errors(request):
    warnings.simplefilter('error')

    request.addfinalizer(lambda *args: warnings.resetwarnings())