"""
Updates all the libraries in ``sunpy/extern``
"""

import os
import json
import shutil
import urllib.request
from parfive import Downloader
from pathlib import Path
from zipfile import ZipFile

SUNPY_EXTERN_DIR = Path(__file__).parent.parent / "sunpy" / "extern"

# "package_name": ["user", "repository", "path_to_file"]
PACKAGES = {
    "appdirs": ["ActiveState", "appdirs", "appdirs.py"],
    "distro": ["python-distro", "distro", "src/distro/distro.py"],
    "inflect": ["jaraco", "inflect", "inflect/__init__.py"],
    "parse": ["r1chardj0n3s", "parse", "parse.py"],
}


def download_package(user: str, repo: str):
    """
    Download the latest version of package using Github release tags.

    Args:
        user: Github user name
        repo: Github repository name
    """
    print(f"Checking {user}/{repo}")
    response = urllib.request.urlopen(f"https://api.github.com/repos/{user}/{repo}")
    if response.status != 200:
        raise ValueError(f"{user}/{repo} does not exist.")

    url = f"https://api.github.com/repos/{user}/{repo}/releases/latest"
    response = urllib.request.urlopen(url)
    if response.status != 200:
        url = f"https://api.github.com/repos/ActiveState/appdirs/tags"
        response = urllib.request.urlopen(url)
        response = json.load(response)
        version = response[0]["name"]
    else:
        response = json.load(response)
        version = response["tag_name"]
    url = f"http://github.com/{user}/{repo}/archive/refs/tags/{version}.zip"

    if not os.path.exists("extern_pkg"):
        os.mkdir("extern_pkg")

    dl = Downloader()
    dl.enqueue_file(url, path="extern_pkg", filename=f"{repo}.zip")
    files = dl.download()
    return f"extern_pkg/{repo}.zip"


def download_github_file(user: str, repo: str, src: Path, dest: Path):
    """
    Download a file from Github.
    Args:
        user: Github user name
        repo: Github repository name
        src: Path to file in the package
        dest: Path where the file should be extracted
    """
    zip_file = download_package(user, repo)
    with ZipFile(zip_file, "r") as f:
        folder = Path(f.namelist()[0]).parts[0]
        ext = f.extract(f"{folder}/{src}", "extern_pkg")
    dest = Path(dest)
    if dest.exists() and dest.is_file():
        os.remove(dest)
    shutil.move(ext, dest)
    os.remove(zip_file)
    shutil.rmtree("extern_pkg")


if __name__ == "__main__":
    for package, (user, repo, src) in PACKAGES.items():
        dest = SUNPY_EXTERN_DIR / f"{package}.py"
        download_github_file(user, repo, src, dest)
