import os.path

testdir = os.path.dirname(os.path.abspath(__file__))

INCLUDE_ONLINE = object()
EXCLUDE_ONLINE = object()
ONLY_ONLINE = object()

try:
    import pytest
except ImportError:
    pytest = None


def main(modulename='', cover=False, show_uncovered_lines=False,
         online=INCLUDE_ONLINE):
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

    cover : bool
        Whether to enable or disable code coverage. The default is False.

    show_uncovered_lines : bool
        Whether to display the line numbers which have not been covered.
        The default is False.

    online : {INCLUDE_ONLINE, EXCLUDE_ONLINE, ONLY_ONLINE}
        A constant from sunpy.tests. EXCLUDE_ONLINE will execute only the
        tests which do not require a working internet connection.
        ONLY_ONLINE will only run the tests which need a connection to the
        internet and INCLUDE_ONLINE, the default, does not restrict the
        number of tests that will be executed in any way.

    """
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
    args = []
    if cover:
        modulepath = os.path.abspath(
            os.path.join(path, os.path.join(os.pardir, os.pardir, modulename)))
        args.extend(['--cov', modulepath])
    if show_uncovered_lines:
        args.extend(['--cov-report', 'term-missing'])
    if online is EXCLUDE_ONLINE:
        args.append('-k-online')
    elif online is ONLY_ONLINE:
        args.extend(['-k', 'online'])
    else:
        if online is not INCLUDE_ONLINE:
            allowed_values = ', '.join(
                'sunpy.tests.' + x for x in [
                    'INCLUDE_ONLINE', 'EXCLUDE_ONLINE', 'ONLY_ONLINE'])
            errmsg = (
                '`online` parameter must have one of the following values: '
                '{0}'.format(allowed_values))
            raise ValueError(errmsg)
    args.append(path)
    return pytest.main(args)
