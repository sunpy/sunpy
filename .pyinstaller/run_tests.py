import os
import sys
import shutil
from pathlib import Path

import pytest

# Skipping these tests that take check the name of the current module
# (which ends up starting with sunpy_tests rather than sunpy).
# asdf/cdf path issue with PyInstaller as well.
SKIP_TESTS = [
    "test_attr_metamagic",
    "test_basic",
    "test_download_cache",
    "test_file_tampered",
    "test_find_dependencies",
    "test_genericmap_basic",
    "test_genericmap_mask",
    "test_hcc_observer_version",
    "test_hgc_100",
    "test_hpc_observer_version",
    "test_main_nonexisting_module",
    "test_main_stdlib_module",
    "test_origin",
    "test_override_file",
    "test_read_cdf",
    "test_same_file_id_different_module",
    "test_saveframe_arr",
    "test_saveframe",
    "test_skip_all",
]
if len(sys.argv) == 3 and sys.argv[1] == "--sunpy-root":
    ROOT = Path(sys.argv[2])
    TEST_ROOT = ROOT / ".pyinstaller" / "sunpy_tests"
    TEST_ROOT.mkdir(parents=True, exist_ok=True)
else:
    # Make sure we don't allow any arguments to be passed - some tests call
    # sys.executable which becomes this script when producing a pyinstaller
    # bundle, but we should just error in this case since this is not the
    # regular Python interpreter.
    if len(sys.argv) > 1 or len(sys.argv) == 1:
        print("Extra arguments passed or missing arguments, exiting early")
        sys.exit(1)
for root, dirnames, files in os.walk(ROOT / "sunpy"):
    for dirname in dirnames:
        final_dir = TEST_ROOT / Path(root).name
        final_dir.mkdir(parents=True, exist_ok=True)
        # We only copy over 'tests' directories, but not sunpy/tests (only
        # sunpy/tests/tests) since that is not just a directory with tests.
        if dirname == "tests" and not root.endswith("sunpy"):
            shutil.copytree(Path(root) / dirname, final_dir, dirs_exist_ok=True)
        else:
            # Create empty __init__.py files so that 'sunpy_tests' still
            # behaves like a single package, otherwise pytest gets confused
            # by the different conftest.py files.
            init_file = final_dir / "__init__.py"
            with open(init_file, "w") as f:
                f.write("#")
    # Copy over all files required to import sunpy that are not in a tests directory.
    for file in files:
        if file in ["conftest.py", "CITATION.rst"]:
            shutil.copy2(Path(root) / file, TEST_ROOT / file)
    # Add the top-level __init__.py file
    with open(TEST_ROOT / "__init__.py", "w") as f:
        f.write("#")
    # Copy the top-level conftest.py
    shutil.copy2(ROOT / "sunpy" / "conftest.py", TEST_ROOT / "conftest.py")
sys.exit(
    pytest.main(
        [
            "sunpy_tests",
            "-s",
            "-vvv",
            "-k " + " and ".join("not " + test for test in SKIP_TESTS),
        ]
    )
)
