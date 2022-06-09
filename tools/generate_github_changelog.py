import subprocess

res = subprocess.run(['pandoc', '--wrap=none', '-t', 'markdown_strict', 'CHANGELOG.rst'], capture_output=True)
print(res.stdout.decode('ascii').split('\n# ')[0])
