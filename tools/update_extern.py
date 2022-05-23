"""
Updates all the libraries in ``sunpy/extern``
"""

import os
import json
import shutil
import tempfile
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
    """
    print(f"Checking {user}/{repo}")
    response = urllib.request.urlopen(f"https://api.github.com/repos/{user}/{repo}")
    if response.status != 200:
        raise ValueError(f"{user}/{repo} does not exist.")

    url = f"https://api.github.com/repos/{user}/{repo}/tags"
    response = urllib.request.urlopen(url)
    if response.status != 200:
        raise ValueError(f"tags for {user}/{repo} does not exist.")
    response = json.load(response)
    version = response[0]["name"]
    url = f"http://github.com/{user}/{repo}/archive/refs/tags/{version}.zip"

    temp_dir = tempfile.mkdtemp()

    print(f"Downloading {user}/{repo}:refs/tags/{version}")
    dl = Downloader()
    dl.enqueue_file(url, path=temp_dir, filename=f"{repo}.zip")
    files = dl.download()
    return files[0]


def download_github_file(user: str, repo: str, src: Path, dest: Path):
    """
    Download a file from Github.
    """
    zip_file = download_package(user, repo)
    temp_dir = tempfile.mkdtemp()
    with ZipFile(zip_file, "r") as f:
        folder = Path(f.namelist()[0]).parts[0]
        ext = f.extract(f"{folder}/{src}", temp_dir)
    dest = Path(dest)
    if dest.exists() and dest.is_file():
        os.remove(dest)
    shutil.move(ext, dest)


if __name__ == "__main__":
    for package, (user, repo, src) in PACKAGES.items():
        download_github_file(user, repo, src, SUNPY_EXTERN_DIR / f"{package}.py")
