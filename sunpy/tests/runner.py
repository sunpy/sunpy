from __future__ import absolute_import, division, print_function

from astropy.tests.runner import TestRunner, keyword


class SunPyTestRunner(TestRunner):
    """
    A runner for SunPy tests, it modifies the arguments to the astropy runner
    to maintain a similar but more SunPy focused CLI.
    """

    # Disable certain astropy flags
    @keyword()
    def remote_data(self, remote_data, kwargs):
        return []

    # Change the docsting on package
    @keyword(priority=10)
    def package(self, package, kwargs):
        """
        package : str, optional
            The name of a specific package to test, e.g. 'map' or 'net.vso'.
            If nothing is specified all default SunPy tests are run.
        """
        return super(SunPyTestRunner, self).package(package, kwargs)

    # Define our own options
    @keyword(False, priority=3)
    def online(self, online, kwargs):
        """
        online : bool, optional
            Enable the online tests if `True` or disable them if `False`.
        """
        if online:
            return ['--remote-data=any']

        return []

    @keyword(False, priority=2)
    def online_only(self, online_only, kwargs):
        """
        online_only: bool, optional
            If `True` only online tests are run.
        """
        r = []
        if online_only:
            r.append('-k remote_data')
            if not kwargs['online']:
                r.append('--remote-data=any')

        return r

    @keyword(False)
    def figure(self, figure, kwargs):
        """
        figure : bool, optional
            Enable the figure tests.
        """
        if not figure and not kwargs['figure_only']:
            return ['-m', 'not figure']

        return []

    @keyword(False)
    def figure_only(self, figure, kwargs):
        """
        figure_only : bool, optional
            Only run the figure tests.
        """
        if figure:
            return ['-m', 'figure']

        return []

    # Define this to change the default value to None
    @keyword()
    def plugins(self, plugins, kwargs):
        """
        plugins : list, optional
            Plugins to be passed to ``pytest.main`` in the ``plugins`` keyword
            argument.
        """
        # Plugins are handled independently by `run_tests` so we define this
        # keyword just for the docstring
        return []
