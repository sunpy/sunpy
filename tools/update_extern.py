"""
Updates all the libraries in ``sunpy/extern``
"""

import os
import shutil
from zipfile import ZipFile
import requests
from pathlib import Path

SUNPY_DIR = Path(__file__).parent.parent


PACKAGES = {
    "appdirs": ["ActiveState", "appdirs", "appdirs.py"],
    "distro": ["python-distro", "distro", "src/distro/distro.py"],
    "inflect": ["jaraco", "inflect", "inflect/__init__.py"],
    "modest_image": ["glue-viz", "glue", "glue/external/modest_image.py"],
    "parse": ["r1chardj0n3s", "parse", "parse.py"],
}


def download_package(package: str):
    """
    Download the latest version of package using Github release tags.

    Args:
        package: The name of the package to download.
    """
    for package in PACKAGES:
        print(f"Checking {PACKAGES[package][0]}/{package}")
        # Get 200 response from github
        response = requests.get(f"https://api.github.com/repos/{PACKAGES[package][0]}/{package}")
        if response.status_code != 200:
            print(f"{PACKAGES[package]}/{package} does not exist.")
            exit()

    author = PACKAGES[package][0]
    url = f"https://api.github.com/repos/{author}/{package}/releases/latest"
    response = requests.get(url)
    version = response.json()["tag_name"]
    url = f"https://github.com/{author}/{package}/archive/refs/tags/{version}.zip"

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


def unzip(folder: str):
    """
    This function unzips the zip file and moves the files to the sunpy/extern folder.

    Args:
        folder: The name of the parent folder of the zip file.
    """
    package = folder.split("-")[0]
    with ZipFile(f"extern_pkg/{package}.zip", "r") as zip_file:
        zip_file.extract(f"{folder}/{PACKAGES[package][1]}", "extern_pkg")
    # os.remove(f"extern_pkg/{package}.zip")


def move(src: Path, dst: Path):
    """
    Move the files from the src to the dst.

    Args:
        src: The path to the files to be moved.
        dst: The path where the files will be moved.
    """
    if os.path.exists(dst):
        os.remove(dst)
    shutil.move(src, dst)


def get_zip_file():
    """
    This function returns a list of parent folders of each zip file inside the temporary folder.
    """

    folders = list()
    for root, dirs, files in os.walk("extern_pkg"):
        for folder in files:
            with ZipFile(f"extern_pkg/{folder}", "r") as zip_file:
                file_name = zip_file.namelist()
                folders.append(file_name[0].split("/")[0])
        folders.sort()
    return folders


if __name__ == "__main__":

    for package in PACKAGES:
        download_package(package)

    folders = get_zip_file()
    for folder in folders:
        unzip(folder)

    for root, dirs, files in os.walk("extern_pkg"):
        for file in files:
            if file.endswith(".py"):
                if file == "__init__.py":
                    package = root.split("/")[-1]
                    move(os.path.join(root, file), os.path.join(
                        SUNPY_DIR, "sunpy", "extern", f"{package}.py"))
                else:
                    move(os.path.join(root, file), os.path.join(SUNPY_DIR, "sunpy", "extern", file))

    shutil.rmtree("extern_pkg")
