from __future__ import absolute_import, division, print_function

import os.path

testdir = os.path.dirname(os.path.abspath(__file__))

try:
    import pytest
except ImportError:
    pytest = None


def main(modulename='', coverage=False, cov_report=False,
         online=False, offline=True, figure=False, verbose=False,
         parallel=0, args=None):
    """
    Execute the test suite of the sunpy package. The parameters may be
    used to restrict the number of tests that will be executed or to
    adjust the output. See the following documentation of the parameters
    for more information.

    Parameters
    ----------
    modulename : str
        The name of the SunPy submodule which will be tested. The default
        is to test the whole sunpy package, i.e. all submodules.

    coverage : bool
        Whether to enable or disable code coverage. The default is False.

    cov_report: string
        Specify if a coverage report should be generated and which one.
        Allowed values: 'html' 'xml' 'annotate' 'term-missing'

    online : bool
        Run the tests that require an internet connection.

    offline: bool
        Run the tests that don't require an internet connection.

    figure: bool
        Include the figure tests in the test run.

    """
    print(modulename)
    if pytest is None:
        raise ImportError("You need to install pytest to run SunPy's tests")

    if not modulename:
        module = __import__('sunpy')
    else:
        module = __import__('sunpy.{0}.tests'.format(modulename), fromlist=[modulename])
    path = None
    for path in module.__path__:
        if os.path.exists(path):
            break
    else:
        raise ImportError(
            'No module named {0!r} in the sunpy package'.format(modulename))
    assert path is not None

    all_args = []
    if coverage:
        print(path, modulename)
        modulepath = os.path.abspath(
            os.path.join(path, os.path.join(os.pardir, os.pardir, modulename)))
        all_args.extend(['--cov', modulepath])
    if cov_report:
        all_args.extend(['--cov-report', cov_report])
    if not online:
        all_args.append('-k-online')
    if not offline:
        all_args.append('-k online')
    if not figure:
        all_args.append('-m not figure')
    all_args.append(path)

    if args:
        all_args.append(args)

    if verbose:
        all_args.append('-v')

    if parallel != 0:
        try:
            import xdist
        except ImportError:
            raise ImportError(
                'Parallel testing requires the pytest-xdist plugin '
                'https://pypi.python.org/pypi/pytest-xdist')

        try:
            parallel = int(parallel)
        except ValueError:
            raise ValueError(
                "parallel must be an int, got {0}".format(parallel))

        all_args.extend(['-n', str(parallel)])

    return pytest.main(all_args)
