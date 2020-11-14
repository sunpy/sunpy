import importlib

from pkg_resources import get_distribution


def find_missing_dependencies(package="sunpy", extras=None):
    """
    List missing dependencies.

    Given a package and, optionally, a tuple of extras, identify any packages
    which should be installed to match the requirements and return any which are
    missing.
    """
    if not extras:
        extras = tuple()

    requires = get_distribution(package).requires(extras=extras)

    missing_requirements = []
    for requirement in requires:
        try:
            get_distribution(requirement)
        except Exception:
            missing_requirements.append(requirement.name)

    return missing_requirements


def missing_dependencies_by_extra(package="sunpy", exclude_extras=None):
    """
    Get all the specified extras for a package and report any missing dependencies.

    This function will also return a "required" item in the dict which is the
    dependencies associated with no extras.
    """
    exclude_extras = exclude_extras or []

    distribution = get_distribution(package)
    extras = distribution.extras

    missing_dependencies = {"required": find_missing_dependencies(package)}
    for extra in extras:
        if extra in exclude_extras:
            continue
        missing_dependencies[extra] = find_missing_dependencies(package, (extra,))

    return missing_dependencies


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


def _self_test_args(*, package=None, online=False, online_only=False, figure_only=False):
    args = ["-W", "ignore"]
    if online:
        args.append("--remote-data=any")
    if online_only:
        args.append("--remote-data=any -m remote_data")
    if package:
        try:
            importlib.import_module(f"sunpy.{package}")
        except ModuleNotFoundError:
            raise ModuleNotFoundError(f"sunpy.{package} was not found.")
        args.extend(["--pyargs", f"sunpy.{package}"])
    else:
        args.extend(["--pyargs", "sunpy"])
    if figure_only:
        args.extend(["-m", "mpl_image_compare"])
    return args


def self_test(*, package=None, online=False, online_only=False, figure_only=False):
    print("\n\n")
    print("Starting sunpy self test...")
    print("Checking for packages needed to run sunpy:")

    missing = missing_dependencies_by_extra(exclude_extras=("asdf", "dask", "dev", "all", "docs"))
    test_missing = missing.pop("tests")
    printed = print_missing_dependencies_report(missing)

    if not printed:
        print("All required and optional sunpy dependencies are installed.")

    if test_missing:
        print("You do not have all the required dependencies installed to run the sunpy test suite.")
        print("NOTE THIS JUST MEANS THE TESTS WILL NOT RUN, NOT THAT SUNPY MISCONFIGURED.")
        print(
            "If you want to run the sunpy tests install the 'tests' extra with `pip install sunpy[all,tests]`")
        return

    import pytest

    print("Starting the sunpy test suite:")
    print()
    args = _self_test_args(package=package, online=online,
                           online_only=online_only, figure_only=figure_only)
    pytest.main(args)
