
def print_missing_dependencies_report(missing, package="sunpy"):
    printed = False
    required_missing = missing.pop("required")
    if required_missing:
        printed = True
        print(f"The following packages are not installed but are required by {package}:")
        for dep in required_missing:
            print(f"* {dep}")
    for extra_name, dependencies in missing.items():
        if not dependencies:
            continue
        printed = True
        print(
            f"The following packages are not installed for the {package}[{extra_name}] requirement:")
        for dep in dependencies:
            print(f"  * {dep}")
    return printed


def _self_test_args(*, package=None, online=False, online_only=False, figure_only=False, verbose=False):
    import importlib

    test_args = ["-W", "ignore"]
    if online:
        test_args.append("--remote-data=any")
    if online_only:
        test_args.append("--remote-data=any -m remote_data")
    if package:
        try:
            importlib.import_module(f"sunpy.{package}")
        except ModuleNotFoundError:
            raise ModuleNotFoundError(f"sunpy.{package} was not found.")
        test_args.extend(["--pyargs", f"sunpy.{package}"])
    else:
        test_args.extend(["--pyargs", "sunpy"])
    if figure_only:
        test_args.extend(["-m", "mpl_image_compare"])
    if verbose:
        test_args.append("-vvv")
    return test_args


def self_test(*, package=None, online=False, online_only=False, figure_only=False, verbose=False):
    from sunpy.util.sysinfo import missing_dependencies_by_extra

    print("Starting sunpy self test...")
    print()
    print("Checking for packages needed to run sunpy:")
    missing = missing_dependencies_by_extra(exclude_extras=["asdf", "dask", "dev", "all", "docs"])
    test_missing = missing.pop("tests")
    printed = print_missing_dependencies_report(missing)
    if not printed:
        print("All required and optional sunpy dependencies are installed.")
    if test_missing:
        print("You do not have all the required dependencies installed to run the sunpy test suite.")
        print()
        print(list(test_missing.keys()))
        print()
        print("If are using conda, you will want to run conda install <package name>")
        print('Otherwise you will want run pip install "sunpy[all,tests]"')
        print()
        return 1
    print()
    print("Starting the sunpy test suite:")
    test_args = _self_test_args(package=package, online=online,
                                online_only=online_only, figure_only=figure_only, verbose=verbose)
    import pytest

    if verbose:
        print(f"Running pytest with arguments: {test_args}")

    return pytest.main(test_args)


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(description="Run the sunpy test suite")
    parser.add_argument("--package", help="The sunpy subpackage to test.")
    parser.add_argument("--online", action="store_true",
                        help="Run the tests that require internet access.")
    parser.add_argument("--online-only", action="store_true",
                        help="Run only the tests that require internet access.")
    parser.add_argument("--figure-only", action="store_true",
                        help="Run only the tests that generate figures.")
    parser.add_argument('--verbose', action='store_true', help="Run in verbose mode")
    test_args = parser.parse_args()
    sys.exit(self_test(package=test_args.package, online=test_args.online,
                       online_only=test_args.online_only, figure_only=test_args.figure_only, verbose=test_args.verbose))
