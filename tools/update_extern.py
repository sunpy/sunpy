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
    "appdirs.py": ["ActiveState", "appdirs", "appdirs.py"],
    "distro.py": ["python-distro", "distro", "src/distro/distro.py"],
    "inflect.py": ["jaraco", "inflect", "inflect/__init__.py"],
    "parse.py": ["r1chardj0n3s", "parse", "parse.py"],
}

LICENSES = {
    "appdirs_license.txt": ["ActiveState", "appdirs", "LICENSE.txt"],
    "distro_license.rst": ["python-distro", "distro", "LICENSE"],
    "inflect_license.txt": ["jaraco", "inflect", "LICENSE"],
    "parse_license.txt": ["r1chardj0n3s", "parse", "LICENSE"],
}


def download_github_file(user: str, repo: str, src: Path, dest: Path):
    """
    Download a file from Github.
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

    url = f"https://raw.githubusercontent.com/{user}/{repo}/{version}/{src}"
    response = urllib.request.urlopen(url)
    if response.status != 200:
        raise ValueError(f"{url} does not exist.")

    with open(dest, "wb") as f:
        print(f"Updating {user}/{repo}:refs/tags/{version}")
        f.write(response.read())
    # zip_file = download_package(user, repo)
    # temp_dir = tempfile.mkdtemp()
    # with ZipFile(zip_file, "r") as f:
    #     folder = Path(f.namelist()[0]).parts[0]
    #     ext = f.extract(f"{folder}/{src}", temp_dir)
    # dest = Path(dest)
    # if dest.exists() and dest.is_file():
    #     os.remove(dest)
    # shutil.move(ext, dest)


if __name__ == "__main__":
    for package, (user, repo, src) in LICENSES.items():
        download_github_file(user, repo, src, SUNPY_EXTERN_DIR / package)
