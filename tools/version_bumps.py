import argparse
import collections
from datetime import datetime, timedelta
from functools import cache
from importlib import metadata

import requests
from packaging.requirements import Requirement
from packaging.version import Version


@cache
def get_package_releases(package):
    """Fetch all versions and their release dates for a package from PyPI"""
    print(f"Querying pypi.org for {package} versions...", end="", flush=True)
    response = requests.get(
        f"https://pypi.org/simple/{package}",
        headers={"Accept": "application/vnd.pypi.simple.v1+json"},
    ).json()
    print("OK")
    file_date = collections.defaultdict(list)
    for f in response["files"]:
        if not f["filename"].endswith(".tar.gz"):
            continue
        if ".dev" in f["filename"] or ".post" in f["filename"]:
            continue
        package_name = response["name"]
        if package_name not in f["filename"]:
            if "-" in package_name:
                # If the package name contains a dash(es), replace it with an underscore
                # to match the filename format used on PyPI
                package_name = package_name.replace("-", "_")
            if "ruamel" in package_name:
                # The filename has a dot, but the package name has an underscore
                package_name = "ruamel.yaml"
            if package_name.capitalize() in f["filename"]:
                package_name = package_name.capitalize()
        try:
            ver = f["filename"].split(package_name + "-")[1].split(".tar.gz")[0]
            version = Version(ver)
        except Exception:
            continue
        release_date = None
        for time_format in ["%Y-%m-%dT%H:%M:%S.%fZ", "%Y-%m-%dT%H:%M:%SZ", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%dT%H:%MZ"]:
            try:
                release_date = datetime.strptime(f["upload-time"], time_format)
                break
            except ValueError:
                continue
        if not release_date:
            continue
        file_date[version].append(release_date)
    release_date = {v: min(file_date[v]) for v in file_date}
    return release_date


def is_version_old(package, version_str, threshold=timedelta(days=365*2)):
    """Check if a specific version of a package is older than the threshold"""
    releases = get_package_releases(package)
    # Find the release date for the specified version
    target_release_date = None
    for v, date in releases.items():
        if v == Version(version_str):
            target_release_date = date
            break
    if not target_release_date:
        print(f"Did not find version {version_str}")
        return True
    now = datetime.now()
    return (now - target_release_date) > threshold


def find_newest_version(package, threshold=timedelta(days=365*2)):
    """Find the oldest available version that is not older than the threshold"""
    releases = get_package_releases(package)
    # Sort releases by date to easily find the most recent ones
    sorted_releases = dict(sorted(releases.items(), key=lambda x: x[0], reverse=False))
    # Find the first release within the threshold
    for v, date in sorted_releases.items():
        if v.is_prerelease or v.micro != 0:
            continue
        now = datetime.now()
        if (now - date) <= threshold:
            return v, date
    return None, None  # No recent versions found


def get_min_version(requirement):
    """
    Extracts the minimum version from a requirement.

    Parameters:
        requirement (packaging.requirements.Requirement): The requirement object.

    Returns:
        str or None: The minimum version if found, otherwise None.
    """
    spec = requirement.specifier
    # Check for specific version constraints
    versions = []
    # Convert the specifier to a string and split by known operators
    spec_str = str(spec)
    operators = ['==', '>=', '>', '<=', '<', '!=']
    # Iterate over each operator to find matching parts in the spec string
    for op in operators:
        if op in spec_str:
            parts = spec_str.split(op)
            if len(parts) > 1 and not any(o in parts[0] for o in operators):
                versions.append((op, parts[1].strip("'\"")))
    # Determine the minimum version based on the extracted parts
    min_version = None
    for op, ver in versions:
        if op == '==' or (op.startswith(('>=', '>')) and not min_version):
            min_version = ver
    return min_version


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
        The values are a nested dictionary with keys being the package names and
        values being the `packaging.requirements.Requirement` objects.
    """
    requirements: list = metadata.requires(package)
    requires_dict = collections.defaultdict(dict)
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
    return requires_dict


def process_dependencies(package, threshold=timedelta(days=365*2)):
    """
    Processes all dependencies to check their versions against the threshold.

    Parameters:
        package (str): The name of the package.
        requirements_dict (dict): Dictionary returned by get_requirements().

    Returns:
        list: A list of dictionaries containing 'package', 'version',
              'group', and 'is_old' for each dependency.
    """
    requirements_dict = get_requirements(package)
    result = collections.defaultdict(list)
    for group, deps in requirements_dict.items():
        for package_name, req in deps.items():
            min_version = get_min_version(req)
            if min_version:
                is_old = is_version_old(package_name, min_version, threshold=threshold)
                result[group].append({
                    'package': package_name,
                    'version': min_version,
                    'group': group,
                    'is_old': is_old
                })
    return result


def output_version_bumps(package, threshold=timedelta(days=365*2)):
    group_deps = process_dependencies(package)
    for group, deps in group_deps.items():
        print(f"\n{group}\n{'-'*len(group)}")
        for dep in deps:
            if dep['is_old']:
                new_version, release_date = find_newest_version(dep['package'], threshold=threshold)
                if new_version is not None:
                    print(f"{dep['package']} should be bumped to {new_version} which was released on {release_date:%Y-%m-%d}")
                else:
                    print(f"Could not find newer version for {dep['package']}")


def main():
    parser = argparse.ArgumentParser(description='Process package dependancies for minimum bumps.')
    # Required positional argument
    parser.add_argument('--packagename', type=str, default="sunpy", help='The name of the package')
    # Optional integer argument with default value
    parser.add_argument('--threshold', type=int, default=720,
                       help='An optional threshold value in days (default: 720)')
    args = parser.parse_args()
    output_version_bumps(args.packagename, threshold=timedelta(days=args.threshold))


if __name__ == "__main__":
    main()
