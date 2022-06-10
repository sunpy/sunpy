"""
Script to render to terminal, the changelog entry for the latest version in markdown so it can be copied and pasted for a GitHub release.

It requires pandoc to be installed on your system.
"""
import subprocess
from pathlib import Path

res = subprocess.run(['pandoc', '--wrap=none', '-t', 'gfm', str(Path(__file__).parent.parent / "CHANGELOG.rst")], capture_output=True)
print(res.stdout.decode('ascii').split('\n# ')[0])
