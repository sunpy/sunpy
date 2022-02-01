"""
This module provides functions to retrieve system information.
"""
import platform
from collections import defaultdict
from importlib.metadata import PackageNotFoundError, version, requires, distribution

from packaging.requirements import Requirement

import sunpy.extern.distro as distro

__all__ = ['system_info', 'find_dependencies', 'missing_dependencies_by_extra']


def get_requirements(package):
    """
    This wraps `importlib.metadata.requires` to not be garbage.

    Parameters
    ----------
    package : str
        Package you want requirements for.

    Returns
    -------
    `dict`
        A dictionary of requirements with keys being the extra requirement group names.
    """
    requirements: list = requires(package)
    requires_dict = defaultdict(list)
    for requirement in requirements:
        req = Requirement(requirement)
        package_name, package_marker = req.name, req.marker
        if package_marker and "extra ==" in str(package_marker):
            group = str(package_marker).split("extra == ")[1].strip('"').strip("'").strip()
            requires_dict[group].append(package_name)
        else:
            requires_dict["required"].append(package_name)
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
        for package in requirements[group]:
            try:
                package_version = version(package)
                installed_requirements[package] = package_version
            except PackageNotFoundError:
                missing_requirements[package] = f"Missing {package}"
    return missing_requirements, installed_requirements


def missing_dependencies_by_extra(package="sunpy", exclude_extras=None):
    """
    Get all the specified extras for a package and report any missing dependencies.

    This function will also return a "required" item in the dict which is the
    dependencies associated with no extras.
    """
    exclude_extras = exclude_extras or []
    requirements = get_requirements(package)
    missing_dependencies = {}
    for group in requirements.keys():
        if group in exclude_extras:
            continue
        missing_dependencies[group] = find_dependencies(package, [group])[0]
    return missing_dependencies


def system_info():
    """
    Prints ones' system info in an "attractive" fashion.
    """
    base_reqs = get_requirements("sunpy")["required"]
    extra_reqs = get_requirements("sunpy")["all"]
    missing_packages, installed_packages = find_dependencies(package="sunpy", extras=["required", "all"])
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
