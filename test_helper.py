# -*- coding: utf-8 -*-
"""
Created on Sat Jun  7 19:36:08 2014

@author: stuart
"""
import os
import sys

sys.path.insert(1, os.path.abspath('astropy_helpers'))
from astropy_helpers.test_helpers import AstropyTest
from astropy_helpers.compat import _fix_user_options


class SunPyTest(AstropyTest):
    description = 'Run the tests for this package'

    user_options = [
        # Package to test
        ('package=', 'P',
         "The name of a specific package to test, e.g. 'io' or 'utils'.  "
         "If nothing is specified, all default tests are run."),
        # Path to test
        #('test-path=', 't',
        # 'Specify a test location by path.  If a relative path to a '
        # '.py file, it is relative to the built package.  If a relative '
        # 'path to a .rst file, it is relative to the docs directory '
        # '(see --docs-path).  May also be an absolute path.'),
        # Print all the things
        ('verbose-results', 'V',
         'Turn on verbose output from pytest.'),
        # plugins to enable
        ('plugins=', 'p',
         'Plugins to enable when running pytest.'),
        # Run only offline tests?
        ('no-online', None,
         'Only run test that do not require a internet connection.'),
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
        #self.test_path = None
        self.verbose_results = False
        self.plugins = None
        self.args = None
        self.no_online = False
        self.coverage = False
        self.cov_report = 'term'
        self.docs_path = None
        self.parallel = 0

    def _validate_required_deps(self):
        pass

    def construct_testing_command(self):
        """
        Build a Python script to run the tests.
        """
        
        cmd_pre = ''  # Commands to run before the test function
        cmd_post = ''  # Commands to run after the test function

        print self.no_online
        print self.coverage
        online = not self.no_online
        print online

        cmd = ('{cmd_pre}{0}; import {1.package_name}, sys; result = ('
               '{1.package_name}.self_test('
               'modulename={1.package!r}, '
               'args={1.args!r}, '
               'verbose={1.verbose_results!r}, '
               'parallel={1.parallel!r}, '
               'online={online!r}, '
               'coverage={1.coverage!r}, '
               'cov_report={1.cov_report!r})); '
               '{cmd_post}'
               'sys.exit(result)')
        return cmd.format('pass', self, online=online, cmd_pre=cmd_pre,
                          cmd_post=cmd_post)
