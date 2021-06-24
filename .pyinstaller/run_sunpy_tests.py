import os
import sys
import shutil

import pytest

# Skipping these tests that take check the name of the current module (which ends up starting with
# sunpy_tests rather than sunpy).
# asdf path issue with PyInstaller as well.

if getattr(sys, 'frozen', False):
    # Running in a bundle
    bundle_dir = sys._MEIPASS

    SKIP_TESTS = [
        'test_saveframe',
        'test_saveframe_arr',
        'test_genericmap_basic',
        'test_genericmap_mask',
        'test_attr_metamagic',
        'test_main_nonexisting_module',
        'test_main_stdlib_module',
        'test_origin',
        'test_find_dependencies',
        'test_missing_dependencies_by_extra',
        'test_hgc_100',
        'test_missing_dependencies_by_extra',
        'test_basic',
        'test_data_manager',
        'test_file_tampered',
        'test_download_cache',
        'test_skip_all',
        'test_override_file',
        'test_same_file_id_different_module']

    sys.exit(pytest.main(['sunpy_tests',
                          '-k ' + ' and '.join('not ' + test for test in SKIP_TESTS)]))
else:
    ROOT = os.path.join(os.path.dirname(__file__), '../')

    for root, dirnames, files in os.walk(os.path.join(ROOT, 'sunpy')):
        for dirname in dirnames:
            final_dir = os.path.relpath(
                os.path.join(
                    root.replace(
                        'sunpy',
                        'sunpy_tests'),
                    dirname),
                ROOT)
            # We only copy over 'tests' directories, but not sunpy/tests (only
            # sunpy/tests/tests) since that is not just a directory with tests.
            if dirname == 'tests' and not root.endswith('sunpy'):
                shutil.copytree(os.path.join(root, dirname), final_dir, dirs_exist_ok=True)
            else:
                # Create empty __init__.py files so that 'sunpy_tests' still
                # behaves like a single package, otherwise pytest gets confused
                # by the different conftest.py files.
                init_filename = os.path.join(final_dir, '__init__.py')
                if not os.path.exists(os.path.join(final_dir, '__init__.py')):
                    os.makedirs(final_dir, exist_ok=True)
                    with open(os.path.join(final_dir, '__init__.py'), 'w') as f:
                        f.write("#")
        # Copy over all conftest.py files
        for file in files:
            if file == 'conftest.py':
                final_file = os.path.relpath(
                    os.path.join(
                        root.replace(
                            'sunpy',
                            'sunpy_tests'),
                        file),
                    ROOT)
                shutil.copy2(os.path.join(root, file), final_file)

        # Add the top-level __init__.py file
        with open(os.path.join('sunpy_tests', '__init__.py'), 'w') as f:
            f.write("#")

        # Copy the top-level conftest.py
        shutil.copy2(os.path.join(ROOT, 'sunpy', 'conftest.py'),
                     os.path.join('sunpy_tests', 'conftest.py'))
