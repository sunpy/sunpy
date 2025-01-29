# Try to use setuptools_scm to get the current version; this is only used
# in development installations from the git repository.
from pathlib import Path

try:
    from setuptools_scm import get_version

    version = get_version(root=Path('../..'), relative_to=__file__)
except ImportError:
    raise ImportError('setuptools_scm not installed')
except Exception as e:
    raise ValueError(f'setuptools_scm broken with {e}')
