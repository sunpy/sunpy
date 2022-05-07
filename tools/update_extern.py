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


def unzip(filename):
    """
    Extract the package.
    """
    with ZipFile(filename, "r") as zip_file:
        zip_file.extractall("tmp")
    # Remove the zip file
    os.remove(filename)


def extract_appdirs(destination):
    """
    Extract the appdirs package.
    """
    # Get appdirs.py from the folder
    print(destination)


if __name__ == "__main__":
    # for package in PACKAGES:
    # # get the url of the package
    # url = get_download_url(package, get_latest_version(package))

    # # download the package
    # download_package(url, package)

    # Open the zip file
    # unzip(f"tmp/{package}.zip")

    # Extract the required files one by one
    extract_appdirs(SUNPY_DIR / "appdirs")
    # shutil.unpack_archive(f"{package}.zip", SUNPY_DIR / "temp")
