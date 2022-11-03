"""
Updates all the libraries in ``sunpy/extern``
"""

import json
import urllib.request
from pathlib import Path

SUNPY_EXTERN_DIR = Path(__file__).parent.parent / "sunpy" / "extern"

# "filename": ["user", "repository", "path_to_file"]
PACKAGES = {
    "appdirs.py": ["ActiveState", "appdirs", "appdirs.py"],
    "appdirs_license.txt": ["ActiveState", "appdirs", "LICENSE.txt"],
    "distro.py": ["python-distro", "distro", "src/distro/distro.py"],
    "distro_license.rst": ["python-distro", "distro", "LICENSE"],
    # Inflect is not here to avoid a dependency on pydantic.
    # Since it only does one specific thing, we should be able to pin the version permanently.
    "parse.py": ["r1chardj0n3s", "parse", "parse.py"],
    "parse_license.txt": ["r1chardj0n3s", "parse", "LICENSE"],
}


def download_github_file(user: str, repo: str, src: Path, dest: Path):
    """
    Download a file from Github.
    """
    try:
        response = urllib.request.urlopen(f"https://api.github.com/repos/{user}/{repo}")
    except urllib.error.HTTPError as e:
        print(f"Error: {e}")
        return

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
        print(f"Updating {dest} to {user}/{repo}:refs/tags/{version}")
        f.write(response.read())


if __name__ == "__main__":
    for package, (user, repo, src) in PACKAGES.items():
        download_github_file(user, repo, src, SUNPY_EXTERN_DIR / package)
