# Try to use setuptools_scm to get the current version; this is only used
# in development installations from the git repository.
from pathlib import Path


"""
Attempt to retrieve the current version of the package using setuptools_scm.
This is only used in development installations from the git repository.

Returns:
    str: The current version of the package.

Raises:
    ImportError: If setuptools_scm is not installed.
    ValueError: If there's an issue with setuptools_scm.
"""
try:
    from setuptools_scm import get_version
    version = get_version(root=Path('../..'), relative_to=__file__)

except ImportError:
    raise ImportError('setuptools_scm is not installed. Please install it to proceed.')

except Exception as e:
    raise ValueError(f'Error while fetching version with setuptools_scm: {e}')

