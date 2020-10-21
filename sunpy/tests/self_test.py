import pytest
from pkg_resources import get_distribution


def find_missing_dependancies(package="sunpy", extras=None):
    """
    List missing dependancies.

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


def missing_dependancies_by_extra(package="sunpy", exclude_extras=None):
    """
    Get all the specified extras for a package and report any missing dependencies.

    This function will also return a "required" item in the dict which is the
    dependencies associated with no extras.
    """
    exclude_extras = exclude_extras or []

    distribution = get_distribution(package)
    extras = distribution.extras

    missing_dependancies = {"required": find_missing_dependancies(package)}
    for extra in extras:
        if extra in exclude_extras:
            continue
        missing_dependancies[extra] = find_missing_dependancies(package, (extra,))

    return missing_dependancies


def print_missing_dependancies_report(missing, package="sunpy"):
    printed = False
    required_missing = missing.pop("required")
    if required_missing:
        printed = True
        print(f"The following packages are not installed but are required by {package}:")
        for dep in required_missing:
            print(f"* {dep}")

    for extra_name, dependancies in missing.items():
        if not dependancies:
            continue

        printed = True
        print(f"The following packages are not installed for the {package}[{extra_name}] requirement:")
        for dep in dependancies:
            print(f"  * {dep}")

    return printed


def self_test(*, online=False):
    print("\n\n")
    print("Starting sunpy self test...")
    print("Checking for installed packages:")

    missing = missing_dependancies_by_extra(exclude_extras=("dev", "all", "docs"))
    test_missing = missing.pop("tests")
    printed = print_missing_dependancies_report(missing)

    if not printed:
        print("All required and optional sunpy dependancies are installed.")

    if test_missing:
        print("You do not have all the required dependancies installed to run the sunpy test suite.")
        print("If you want to run the sunpy tests install the 'tests' extra with `pip install sunpy[tests]`")
        return

    print("Starting the sunpy test suite:")
    print()
    args = ["--pyargs", "sunpy", "-W", "ignore"]
    if online:
        args.append("--remote-data=any")

    pytest.main(args)
