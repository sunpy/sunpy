"""
This module provides functions to retrieve system information.
"""
import platform
from collections import defaultdict
from importlib.metadata import PackageNotFoundError, version, requires, distribution

from packaging.requirements import Requirement

import sunpy.extern.distro as distro

__all__ = ['system_info', 'find_dependencies', 'missing_dependencies_by_extra']


def get_requirements(package, names_only=False):
    """
    This wraps `importlib.metadata.requires` to not be garbage.

    Parameters
    ----------
    package : str
        Package you want requirements for.
    names_only : boolean
        Whether to return just the requirement package names or more details, including the extra group

    Returns
    -------
    `dict`
        A dictionary of requirements with keys being the extra requirement group names.
        The values are either requirement package names, or a nested dictionary with keys being the
        package names and values being the `packaging.requirements.Requirement` objects.
    """
    requirements: list = requires(package)
    if names_only:
        requires_dict = defaultdict(list)
    else:
        requires_dict = defaultdict(dict)
    for requirement in requirements:
        req = Requirement(requirement)
        package_name, package_marker = req.name, req.marker
        if package_marker and "extra ==" in str(package_marker):
            group = str(package_marker).split("extra == ")[1].strip('"').strip("'").strip()
        else:
            group = "required"
        # De-duplicate (the same package could appear more than once in the extra == 'all' group)
        if package_name in requires_dict[group]:
            continue
        if names_only:
            requires_dict[group].append(package_name)
        else:
            requires_dict[group][package_name] = req
    return requires_dict


def find_dependencies(package="sunpy", extras=None):
    """
    List installed and missing dependencies.

    Given a package and, optionally, a tuple of extras, identify any packages
    which should be installed to match the requirements and return any which are
    missing.
    """
    requirements = get_requirements(package)
    installed_requirements = {}
    missing_requirements = {}
    extras = extras or ["required"]
    for group in requirements:
        if group not in extras:
            continue
        for package, package_details in requirements[group].items():
            try:
                package_version = version(package)
                installed_requirements[package] = package_version
            except PackageNotFoundError:
                missing_requirements[package] = f"Missing {package_details}"
    return missing_requirements, installed_requirements


def missing_dependencies_by_extra(package="sunpy", exclude_extras=None):
    """
    Get all the specified extras for a package and report any missing dependencies.

    This function will also return a "required" item in the dict which is the
    dependencies associated with no extras.
    """
    exclude_extras = exclude_extras or []
    requirements = get_requirements(package, True)
    missing_dependencies = {}
    for group in requirements.keys():
        if group in exclude_extras:
            continue
        missing_dependencies[group] = find_dependencies(package, [group])[0]
    return missing_dependencies


def get_extra_groups(groups, exclude_extras):
    return list(set(groups) - set(exclude_extras))


def system_info():
    """
    Prints ones' system info in an "attractive" fashion.
    """
    requirements = get_requirements("sunpy", True)
    groups = [*requirements.keys()]
    extra_groups = get_extra_groups(groups, ['all', 'dev'])
    base_reqs = requirements['required']
    extra_reqs = requirements['all']
    missing_packages, installed_packages = find_dependencies(package="sunpy", extras=extra_groups)
    extra_prop = {"System": platform.system(),
                  "Arch": f"{platform.architecture()[0]}, ({platform.processor()})",
                  "Python": platform.python_version(),
                  "sunpy": version("sunpy")}
    sys_prop = {**installed_packages, **missing_packages, **extra_prop}
    print("==============================")
    print("sunpy Installation Information")
    print("==============================")
    print()
    print("General")
    print("#######")
    if sys_prop['System'] == "Linux":
        print(f"OS: {distro.name()} ({distro.version()}, Linux {platform.release()})")
    elif sys_prop['System'] == "Darwin":
        print(f"OS: Mac OS {platform.mac_ver()[0]}")
    elif sys_prop['System'] == "Windows":
        print(f"OS: Windows {platform.release()} {platform.version()}")
    else:
        print("Unknown OS")
    for sys_info in ['Arch', 'sunpy']:
        print(f'{sys_info}: {sys_prop[sys_info]}')
    print(f'Installation path: {distribution("sunpy")._path}')
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
