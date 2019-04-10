import os

from astropy.tests.runner import TestRunner, keyword


class SunPyTestRunner(TestRunner):
    """
    A runner for SunPy tests, it modifies the arguments to the astropy runner
    to maintain a similar but more SunPy focused CLI.
    """
    # Disable certain astropy flags
    @keyword()
    def remote_data(self, remote_data, kwargs):
        return NotImplemented  # pragma: no cover

    # Change the docsting on package
    @keyword(priority=10)
    def package(self, package, kwargs):
        """
        package : `str`, optional
            The name of a specific package to test, e.g. 'map' or 'net.vso'.
            If nothing is specified all default SunPy tests are run.
        """
        return super().package(package, kwargs)

    # Define our own options
    @keyword(False, priority=3)
    def online(self, online, kwargs):
        """
        online : `bool`, optional
            Enable the online tests if `True` or disable them if `False`.
        """
        if online:
            return ['--remote-data=any']

        return []

    @keyword(False, priority=2)
    def online_only(self, online_only, kwargs):
        """
        online_only: `bool`, optional
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
        figure : `bool`, optional
            Enable the figure tests.
        """
        if not figure and not kwargs['figure_only']:
            return ['-m', 'not figure']

        return []

    @keyword(False)
    def figure_only(self, figure, kwargs):
        """
        figure_only : `bool`, optional
            Only run the figure tests.
        """
        if figure:
            return ['-m', 'figure']

        return []

    @keyword()
    def figure_dir(self, figure_dir, kwargs):
        """
        figure_tests : `str`, optional
            Set the output directory for figure test images and hashes.
        """
        # If our test path is outside of our base dir (docs) then we have to
        # skip sending --figure_dir as we will not have hit the conftest.py
        # file.
        if kwargs['test_path'] and self.base_path not in kwargs['test_path']:
            return []

        if figure_dir:
            return ['--figure_dir', figure_dir]

        return []

    # Define this to change the default value to None
    @keyword()
    def plugins(self, plugins, kwargs):
        """
        plugins : `list`, optional
            Plugins to be passed to ``pytest.main`` in the ``plugins`` keyword
            argument.
        """
        # Plugins are handled independently by `run_tests` so we define this
        # keyword just for the docstring
        return []

    @keyword()
    def coverage(self, coverage, kwargs):
        if coverage:
            coveragerc = os.path.join(self.base_path, "tests", "coveragerc")
            ret = []
            for path in self.package_path:
                ret += ["--cov", path, "--cov-config", coveragerc]
            return ret

        return []

    @keyword()
    def cov_report(self, cov_report, kwargs):
        if kwargs['coverage'] and cov_report:
            a = [cov_report] if isinstance(cov_report, str) else []
            return ['--cov-report'] + a

        return []

    @keyword()
    def docs_path(self, docs_path, kwargs):
        """
        Ignore the docs directory if we have set test_path.
        """
        if kwargs['test_path']:
            return []
        else:
            return super().docs_path(docs_path, kwargs)
