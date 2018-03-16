# -*- coding: utf-8 -*-
"""
Created on Sat Jun  7 19:36:08 2014

@author: Stuart Mumford

This file is designed to be imported and ran only via setup.py, hence it's
dependency on astropy_helpers which will be available in that context.
"""
from __future__ import absolute_import, division, print_function

import os

from astropy_helpers.commands.test import AstropyTest
from astropy_helpers.compat import _fix_user_options


class SunPyTest(AstropyTest):
    description = 'Run the tests for this package'

    user_options = [
        # Package to test
        ('package=', 'P',
         "The name of a specific package to test, e.g. 'io' or 'utils'.  "
         "If nothing is specified, all default tests are run."),
        # Print all the things
        ('verbose-results', 'V',
         'Turn on verbose output from pytest.'),
        # plugins to enable
        ('plugins=', 'p',
         'Plugins to enable when running pytest.'),
        # Run online tests?
        ('online', 'R',
         'Also run tests that do require a internet connection.'),
        # Run only online tests?
        ('online-only', None,
         'Only run test that do require a internet connection.'),
        # Run tests that check figure generation
        ('figure', None,
         'Run tests that compare figures against stored hashes.'),
        # Run only tests that check figure generation
        ('figure-only', None,
         'Only run tests that compare figures against stored hashes.'),
        # Calculate test coverage
        ('coverage', 'c',
         'Create a coverage report. Requires the coverage package.'),
        ('cov-report=', None,
         'Specify the type of coverage report to generate. (Default terminal)'),
        # Run tests in parallel
        ('parallel=', 'j',
         'Run the tests in parallel on the specified number of '
         'CPUs.  If negative, all the cores on the machine will be '
         'used.  Requires the pytest-xdist plugin.'),
        # Pass additional cli args to pytest
        ('args=', 'a',
         'Additional arguments to be passed to pytest.')
    ]

    user_options = _fix_user_options(user_options)
    package_name = ''

    def initialize_options(self):
        self.package = ''
        self.verbose_results = False
        self.plugins = None
        self.args = None
        self.online = False
        self.online_only = False
        self.figure = False
        self.figure_only = False
        self.coverage = False
        self.cov_report = 'term' if self.coverage else None
        self.docs_path = os.path.abspath('docs')
        self.parallel = 0
        self.temp_root = None

    def _validate_required_deps(self):
        """
        This method checks that any required modules are installed before
        running the tests.
        """
        try:
            import sunpy
        except ImportError:
            raise ImportError(
                "The 'test' command requires the sunpy package to be "
                "installed and importable.")

    def generate_testing_command(self):
        """
        Build a Python script to run the tests.
        """

        self._validate_required_deps()  # checks if modules are installed

        cmd_pre = ''  # Commands to run before the test function
        cmd_post = ''  # Commands to run after the test function

        if self.coverage:
            pre, post = self._generate_coverage_commands()
            cmd_pre += pre
            cmd_post += post

        online = self.online
        offline = not self.online_only
        figure = self.figure
        figure_only = self.figure_only

        cmd = ('{cmd_pre}{0}; import {1.package_name}, sys; result = ('
               '{1.package_name}.self_test('
               'modulename={1.package!r}, '
               'args={1.args!r}, '
               'verbose={1.verbose_results!r}, '
               'parallel={1.parallel!r}, '
               'online={online!r}, '
               'offline={offline!r}, '
               'figure={figure!r}, '
               'figure_only={figure_only!r}, '
               'cov_report={1.cov_report!r})); '
               '{cmd_post}'
               'sys.exit(result)')
        x = cmd.format('pass',
                       self,
                       online=online,
                       offline=offline,
                       figure=figure,
                       figure_only=figure_only,
                       cmd_pre=cmd_pre,
                       cmd_post=cmd_post)
        return x
