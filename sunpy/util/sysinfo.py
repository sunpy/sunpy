"""
This module provides functions to retrieve system information.
"""
import platform

from pkg_resources import get_distribution

from sunpy.extern.distro import linux_distribution

__all__ = ['system_info', 'find_dependencies', 'missing_dependencies_by_extra']


def find_dependencies(package="sunpy", extras=None):
    """
    List installed and missing dependencies.

    Given a package and, optionally, a tuple of extras, identify any packages
    which should be installed to match the requirements and return any which are
    missing.
    """
    if not extras:
        extras = tuple()
    requires = get_distribution(package).requires(extras=extras)
    installed_requirements = {}
    missing_requirements = {}
    for requirement in requires:
        try:
            package = get_distribution(requirement)
            installed_requirements[package.project_name.lower()] = package.version
        except Exception:
            missing_requirements[requirement.name.lower()] = f"Missing, need {requirement}"
    return missing_requirements, installed_requirements


def missing_dependencies_by_extra(package="sunpy", exclude_extras=None):
    """
    Get all the specified extras for a package and report any missing dependencies.

    This function will also return a "required" item in the dict which is the
    dependencies associated with no extras.
    """
    exclude_extras = exclude_extras or []
    distribution = get_distribution(package)
    extras = distribution.extras
    missing_dependencies = {"required": find_dependencies(package)[0]}
    for extra in extras:
        if extra in exclude_extras:
            continue
        missing_dependencies[extra] = find_dependencies(package, (extra,))[0]
    return missing_dependencies


def system_info():
    """
    Prints ones' system info in an "attractive" fashion.
    """
    base_reqs = get_distribution("sunpy").requires()
    base_reqs = {base_req.name.lower() for base_req in base_reqs}
    extra_reqs = get_distribution("sunpy").requires(extras=["all"])
    extra_reqs = sorted({extra_req.name.lower() for extra_req in extra_reqs}.difference(base_reqs))

    missing_packages, installed_packages = find_dependencies(package="sunpy", extras=["all"])
    extra_prop = {"System": platform.system(),
                  "Arch": f"{platform.architecture()[0]}, ({platform.processor()})",
                  "Python": platform.python_version(),
                  "sunpy": get_distribution("sunpy").version}
    sys_prop = {**installed_packages, **missing_packages, **extra_prop}

    print("==============================")
    print("sunpy Installation Information")
    print("==============================")
    print()
    print("General")
    print("#######")
    if sys_prop['System'] == "Linux":
        distro = " ".join(linux_distribution())
        print(f"OS: {distro} (Linux {platform.release()})")
    elif sys_prop['System'] == "Darwin":
        print(f"OS: Mac OS {platform.mac_ver()[0]}")
    elif sys_prop['System'] == "Windows":
        print(f"OS: Windows {platform.release()} {platform.version()}")
    else:
        print("Unknown OS")
    for sys_info in ['Arch', 'sunpy']:
        print(f'{sys_info}: {sys_prop[sys_info]}')
    print(f'Installation path: {get_distribution("sunpy").location}')
    print()
    print("Required Dependencies")
    print("#####################")
    for req in base_reqs:
        print(f'{req}: {sys_prop[req]}')
    print()
    print("Optional Dependencies")
    print("#####################")
    for extra_req in extra_reqs:
        print(f'{extra_req}: {sys_prop[extra_req]}')
