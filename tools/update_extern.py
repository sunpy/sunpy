"""
Updates all the libraries in ``sunpy/extern``
"""

import os
import shutil
import requests
from tqdm import tqdm
from pathlib import Path
from zipfile import ZipFile

SUNPY_EXTERN_DIR = Path(__file__).parent / "extern"

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
        package: The name of the package to download.
    """
    print(f"Checking {user}/{repo}")
    response = requests.get(f"https://api.github.com/repos/{user}/{repo}")
    if response.status_code != 200:
        print(f"{user}/{repo} does not exist.")
        exit()

    url = f"https://api.github.com/repos/{user}/{repo}/releases/latest"
    response = requests.get(url)
    if response.status_code != 200:
        url = f"https://api.github.com/repos/{user}/{repo}/tags"
        response = requests.get(url)
        version = response.json()[0]["name"]
    else:
        version = response.json()["tag_name"]
    url = f"https://github.com/{user}/{repo}/archive/refs/tags/{version}.zip"

    if not os.path.exists("extern_pkg"):
        os.mkdir("extern_pkg")

    response = requests.get(url, stream=True)
    with open(f"extern_pkg/{repo}.zip", "wb") as f:
        print(f"Downloading {repo}")
        for chunk in tqdm(response.iter_content(chunk_size=1024)):
            if chunk:
                f.write(chunk)
            f.flush()
    return f"extern_pkg/{repo}.zip"


def move(src: Path, dest: Path):
    """
    Move the files from the src to the dst.

    Args:
        src: The path to the files to be moved.
        dest: The path where the files will be moved.
    """
    src = Path(src)
    dest = Path(dest)
    if dest.exists():
        os.remove(dest)
    shutil.move(src, dest)


def get_zip_file():
    """
    This function returns name of the parent folders of zip file inside the temporary folder.
    """

    for root, dirs, files in os.walk("extern_pkg"):
        for folder in files:
            with ZipFile(f"extern_pkg/{folder}", "r") as zip_file:
                file_name = zip_file.namelist()
                return (file_name[0].split("/")[0])


def download_github_file(user: str, repo: str, src: Path, dest: Path):
    zip_file = download_package(user, repo)
    with ZipFile(zip_file, "r") as f:
        folder = Path(f.namelist()[0]).parts[0]
        ext = f.extract(Path(folder) / src, "extern_pkg")
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
