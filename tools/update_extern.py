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
    if not os.path.exists("tmp"):
        os.mkdir("tmp")

    response = requests.get(url, stream=True)
    with open(f"tmp/{package}.zip", "wb") as f:
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
        zip_file.extractall("tmp")
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
    # shutil.move(src, dst)
    print(f"Moving {src} to {dst}")


def extract_appdirs(destination):
    """
    Extract the appdirs package.
    """
    # Print the path of appdirs.py

    if os.path.exists(f"{destination}/appdirs.py"):
        print(f"{destination}/appdirs.py")
        # Move the file to the sunpy/extern directory
        move_files(f"{destination}/appdirs.py", SUNPY_DIR / "extern" / "appdirs.py")


def extract_distro(destination):
    """
    Extract the distro package.
    """
    # Print the path of distro.py
    if os.path.exists(f"{destination}/src/distro/distro.py"):
        print(f"{destination}/src/distro/distro.py")


def extract_inflect(destination):
    """
    Extract the inflect package.
    """
    if os.path.exists(f"{destination}/inflect/__init__.py"):
        print(f"{destination}/inflect/__init__.py")


def extract_parse(destination):
    """
    Extract the parse package.
    """
    if os.path.exists(f"{destination}/parse.py"):
        print(f"{destination}/parse.py")


if __name__ == "__main__":
    # for package in PACKAGES:
    # # get the url of the package
    # url = get_download_url(package, get_latest_version(package))

    # # download the package
    # download_package(url, package)

    # Open the zip file
    # unzip(f"tmp/{package}.zip")

    # Extract the required files one by one

    folder_name = list()
    for root, dirs, files in os.walk("tmp"):
        folder_name = dirs
        break
    folder_name.sort()

    # Extract the appdirs package
    extract_appdirs(f"tmp/{folder_name[0]}")

    # Extract the distro package
    extract_distro(f"tmp/{folder_name[1]}")

    # Extract the inflect package
    extract_inflect(f"tmp/{folder_name[2]}")

    # Extract the parse package
    extract_parse(f"tmp/{folder_name[3]}")
