#!/bin/env python

"""
This script is used to update the files in `sunpy/extern` while making new releases.
"""

import os
import shutil
from zipfile import ZipFile
import requests
from pathlib import Path

# get parent directory of sunpy
SUNPY_DIR = Path(__file__).parent.parent


# List of packages to update
PACKAGES = ["appdirs", "distro", "inflect", "parse"]

# Dictionary to store the authors of the packages
AUTHORS = {
    "appdirs": "ActiveState",
    "distro": "python-distro",
    "inflect": "jaraco",
    "parse": "r1chardj0n3s",
}


def get_latest_version(package):
    """
    Get the latest version of the package.
    """
    url = f"https://api.github.com/repos/{AUTHORS[package]}/{package}/releases/latest"
    response = requests.get(url)
    return response.json()["tag_name"]


def get_download_url(package, version):
    """
    Get the download url for the package.
    """
    url = f"https://github.com/{AUTHORS[package]}/{package}/archive/refs/tags/{version}.zip"
    return url


def download_package(url, package):
    """
    Download the package.
    """
    # Make a temporaty directory
    if not os.path.exists("extern_pkg"):
        os.mkdir("extern_pkg")

    response = requests.get(url, stream=True)
    with open(f"extern_pkg/{package}.zip", "wb") as f:
        print(f"Downloading {package}")
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
            f.flush()
        f.close()


def unzip(filename):
    """
    Extract the package.
    """
    with ZipFile(filename, "r") as zip_file:
        zip_file.extractall("extern_pkg")
    # Remove the zip file
    os.remove(filename)


def move_files(src: Path, dst: Path):
    """
    Move the files from the src to the dst.
    """
    # if dst already exists, remove it
    if os.path.exists(dst):
        print(f"Removing {dst}")
    # move the files
    shutil.move(src, dst)
    print(f"Moving {src} to {dst}")


def update_extern(destination):
    """
    Update the files in the sunpy/extern directory.
    """
    # If the destination contains the appdirs package, extract the appdirs.py file
    if destination.split("-")[0] == "appdirs":
        if os.path.exists(f"extern_pkg/{destination}/appdirs.py"):
            # Move the file to the sunpy/extern directory
            move_files(f"extern_pkg/{destination}/appdirs.py", SUNPY_DIR / "sunpy" / "extern" / "appdirs.py")

    # If the destination contains the distro package, extract the distro.py file
    if destination.split("-")[0] == "distro":
        # Print the path of distro.py
        if os.path.exists(f"extern_pkg/{destination}/src/distro/distro.py"):
            # Move the file to the sunpy/extern directory
            move_files(f"extern_pkg/{destination}/src/distro/distro.py",
                       SUNPY_DIR / "sunpy" / "extern" / "distro.py")

    # If the destination contains the inflect package, extract the inflect.py file
    if destination.split("-")[0] == "inflect":
        if os.path.exists(f"extern_pkg/{destination}/inflect/__init__.py"):
            # Rename the file to inflect.py
            os.rename(f"extern_pkg/{destination}/inflect/__init__.py",
                      f"extern_pkg/{destination}/inflect/inflect.py")
            # Move the file to the sunpy/extern directory
            move_files(f"extern_pkg/{destination}/inflect/inflect.py",
                       SUNPY_DIR / "sunpy" / "extern" / "inflect.py")

    # If the destination contains the parse package, extract the parse.py file
    if destination.split("-")[0] == "parse":
        if os.path.exists(f"extern_pkg/{destination}/parse.py"):
            # Move the file to the sunpy/extern directory
            move_files(f"extern_pkg/{destination}/parse.py", SUNPY_DIR / "sunpy" / "extern" / "parse.py")

    # Remove the directory
    shutil.rmtree(f"extern_pkg/{destination}")


if __name__ == "__main__":
    for package in PACKAGES:
        # get the url of the package
        url = get_download_url(package, get_latest_version(package))

        # download the package
        download_package(url, package)

        # Open the zip file
        unzip(f"extern_pkg/{package}.zip")

    folders = list()
    for root, dirs, files in os.walk("extern_pkg"):
        folders = dirs
        break

    # Sort the packages in alphabetical order
    folders.sort()

    # Extract the files
    for folder in folders:
        update_extern(folder)

    # Remove the temporary directory
    shutil.rmtree("extern_pkg")
