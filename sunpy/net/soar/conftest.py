import os

import pytest


@pytest.fixture(scope="session", autouse=True)
def _hide_parfive_progress(request):
    """
    Set the PARFIVE_HIDE_PROGRESS to hide the parfive progress bar in tests.
    """
    os.environ["PARFIVE_HIDE_PROGRESS"] = "True"
    yield
    del os.environ["PARFIVE_HIDE_PROGRESS"]
