"""
This file is designed to be imported and ran only via "setup.py", hence it's
dependency on astropy_helpers which will be available in that context.
"""
import os
import copy

from astropy_helpers.commands.test import AstropyTest


class SunPyTest(AstropyTest):
    description = 'Run the tests for this package'

    user_options = copy.copy(AstropyTest.user_options)

    user_options.remove(('remote-data=', 'R',
                         'Run tests that download remote data. Should be '
                         'one of none/astropy/any (defaults to none).'))

    user_options += [('online', 'R',
                      'Also run tests that do require a internet connection.'),
                     ('online-only', None,
                      'Only run test that do require a internet connection.'),
                     ('cov-report=', None,
                      'How to display the coverage report, should be either "html" or "term"'),
                     ('figure', None,
                      'Run the figure tests.'),
                     # Run only tests that check figure generation
                     ('figure-only', None,
                      'Only run tests that compare figures against stored hashes.')]

    package_name = ''

    def initialize_options(self):
        super().initialize_options()
        self.online = False
        self.online_only = False
        self.figure = False
        self.figure_only = False
        self.cov_report = True

    def _generate_coverage_commands(self):
        cmd_pre = ''  # Commands to run before the test function

        # patch the .coverage file so the paths are correct to the directory
        # setup.py was run in rather than the temporary directory.
        cwd = os.path.abspath(".")
        cmd_post = ('from sunpy.tests.helpers import _patch_coverage; '
                    'import os; '
                    'test_dir = os.path.abspath("."); '
                    f'_patch_coverage(test_dir, "{cwd}"); ')

        # Make html report the default and make pytest-cov save it to the
        # source directory not the temporary directory.
        if self.cov_report and (isinstance(self.cov_report, bool) or "html" in self.cov_report):
            html_cov = os.path.join(os.path.abspath("."), "htmlcov")
            self.cov_report = f'html:{html_cov}'
        else:
            self.cov_report = self.cov_report

        return cmd_pre, cmd_post

    def generate_testing_command(self):
        """
        Build a Python script to run the tests.
        """

        cmd_pre = ''  # Commands to run before the test function
        cmd_post = ''  # Commands to run after the test function

        if self.coverage:
            pre, post = self._generate_coverage_commands()
            cmd_pre += pre
            cmd_post += post

        cmd = ('{cmd_pre}{0}; import {1.package_name}, sys; result = ('
               '{1.package_name}.self_test('
               'package={1.package!r}, '
               'test_path={1.test_path!r}, '
               'args={1.args!r}, '
               'coverage={1.coverage!r}, '
               'cov_report={1.cov_report!r}, '
               'plugins={1.plugins!r}, '
               'verbose={1.verbose_results!r}, '
               'pastebin={1.pastebin!r}, '
               'online={1.online!r}, '
               'online_only={1.online_only!r}, '
               'figure={1.figure!r}, '
               'figure_only={1.figure_only!r}, '
               'figure_dir="{figure_dir}", '
               'pep8={1.pep8!r}, '
               'pdb={1.pdb!r}, '
               'open_files={1.open_files!r}, '
               'parallel={1.parallel!r}, '
               'docs_path={1.docs_path!r}, '
               'skip_docs={1.skip_docs!r}, '
               'repeat={1.repeat!r})); '
               '{cmd_post}'
               'sys.exit(result)')
        return cmd.format('pass',
                          self,
                          figure_dir=os.path.join(os.path.abspath('.'), "figure_test_images"),
                          cmd_pre=cmd_pre,
                          cmd_post=cmd_post)
