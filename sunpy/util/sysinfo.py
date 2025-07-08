"""
This module provides functions to retrieve system information.
"""
import sys
import platform
from collections import defaultdict
from importlib.metadata import PackageNotFoundError, version, requires, distribution

from packaging.markers import Marker
from packaging.requirements import Requirement

import sunpy.extern.distro as distro
from sunpy.util import warn_user

__all__ = ['system_info', 'find_dependencies', 'missing_dependencies_by_extra']


def _trusted_version(package_name):
    """
    If the package has already been imported, trust its __version__ attribute
    over the version retrievable by importlib.
    """
    if package_name in sys.modules:
        package = sys.modules[package_name]
        if hasattr(package, "__version__"):
            return package.__version__
    return version(package_name)


def get_requirements(package, *, expand_groups=False):
    """
    This wraps `importlib.metadata.requires` to not be garbage.

    Parameters
    ----------
    package : str
        Package you want requirements for.

    expand_groups : bool
        If `True`, expand any requirement that is a group of the specified package.

    Returns
    -------
    `dict`
        A dictionary of requirements with keys being the extra requirement group names.
        The values are a nested dictionary with keys being the package names and
        values being the `packaging.requirements.Requirement` objects.
    """
    requirements: list = requires(package)
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
        requires_dict[group][package_name] = req

    if expand_groups:
        # First expand each self-referential package requirement into the individual groups
        for group, group_reqs in requires_dict.items():
            if package in group_reqs and (extras := group_reqs[package].extras):
                requires_dict[group].update(dict.fromkeys(extras))
                del group_reqs[package]

        # Resolve each group, recursing as necessary
        for group in requires_dict:
            _resolve_group(group, requires_dict)

    return requires_dict


def _resolve_group(group, requires_dict):
    """
    Return a fully resolved list of requirements for a group.
    If the group requires another group, recurse to fully resolve that group first.
    The dictionary is permanently updated with the fully resolved requirements.
    """
    for req_name, req in requires_dict[group].copy().items():
        if req is None:  # req==None means that req_name is an unresolved group
            requires_dict[group].update(_resolve_group(req_name, requires_dict))
            del requires_dict[group][req_name]
    return requires_dict[group]


def resolve_requirement_versions(package_versions):
    """
    Resolves a list of requirements for the same package.

    Given a list of package details in the form of `packaging.requirements.Requirement`
    objects, combine the specifier, extras, url and marker information to create
    a new requirement object.
    """
    resolved = Requirement(str(package_versions[0]))

    for package_version in package_versions[1:]:
        resolved.specifier = resolved.specifier & package_version.specifier
        resolved.extras = resolved.extras.union(package_version.extras)
        resolved.url = resolved.url or package_version.url
        if resolved.marker and package_version.marker:
            resolved.marker = Marker(f"{resolved.marker} or {package_version.marker}")
        elif package_version.marker:
            resolved.marker = package_version.marker

    return resolved


def format_requirement_string(requirement):
    formatted_string = f"Missing {requirement}"
    formatted_string = formatted_string.replace("or extra ==", "or").strip()
    return formatted_string


def find_dependencies(package="sunpy", extras=None):
    """
    List installed and missing dependencies.

    Given a package and, optionally, a tuple of extras, identify any packages
    which should be installed to match the requirements and return any which are
    missing.
    """
    requirements = get_requirements(package, expand_groups=True)
    installed_requirements = {}
    missing_requirements = defaultdict(list)
    extras = extras or ["required"]
    for group in requirements:
        if group not in extras:
            continue
        for package, package_details in requirements[group].items():
            try:
                package_version = _trusted_version(package)
                installed_requirements[package] = package_version
            except PackageNotFoundError:
                missing_requirements[package].append(package_details)
    for package, package_versions in missing_requirements.items():
        missing_requirements[package] = format_requirement_string(
            resolve_requirement_versions(package_versions))
    return missing_requirements, installed_requirements


def _warn_missing_deps(extras):
    """
    Warn a user if they are missing dependencies defined in a given extras.
    """
    if (deps := find_dependencies(package="sunpy", extras=extras)):
        missing_deps = [deps[0][key].split(";")[0].removeprefix("Missing ") for key in deps[0].keys()]
        if missing_deps:
            warn_user(f"Importing sunpy.{extras} without its extra dependencies may result in errors.\n"
                      f"The following packages are not installed:\n{missing_deps}\n"
                      f"To install sunpy with these dependencies use `pip install sunpy[{extras}]` "
                      f"or `pip install sunpy[all]` for all extras. \n"
                      "If you installed sunpy via conda, please report this "
                      "to the community channel: https://matrix.to/#/#sunpy:openastronomy.org"
                      )


def missing_dependencies_by_extra(package="sunpy", exclude_extras=None):
    """
    Get all the specified extras for a package and report any missing dependencies.

    This function will also return a "required" item in the dict which is the
    dependencies associated with no extras.
    """
    exclude_extras = exclude_extras or []
    requirements = get_requirements(package)  # groups do not need to be expanded here
    missing_dependencies = {}
    for group in requirements.keys():
        if group in exclude_extras:
            continue
        missing_dependencies[group] = find_dependencies(package, [group])[0]
    return missing_dependencies


def get_extra_groups(groups, exclude_extras):
    return list(set(groups) - set(exclude_extras))


def get_keys_list(dictionary, sort=True):
    keys = [*dictionary.keys()]
    if sort:
        return sorted(keys)
    return keys


def system_info():
    """
    Prints ones' system info in an "attractive" fashion.
    """
    requirements = get_requirements("sunpy", expand_groups=True)
    groups = get_keys_list(requirements)
    extra_groups = get_extra_groups(groups, ['all', 'dev'])
    base_reqs = get_keys_list(requirements['required'])
    extra_reqs = get_keys_list(requirements['all'])
    missing_packages, installed_packages = find_dependencies(package="sunpy", extras=extra_groups)
    extra_prop = {"System": platform.system(),
                  "Arch": f"{platform.architecture()[0]}, ({platform.processor()})",
                  "Python": platform.python_version(),
                  "sunpy": _trusted_version("sunpy")}
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
