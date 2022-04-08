import os
import re
import glob
import subprocess


def inc(obj, key):
    """Increment a dictionary[key] value."""
    try:
        obj[key] += 1
    except KeyError:
        obj[key] = 1


# Load PR number, comment and bot name
try:
    pr = os.environ['PR_NUMBER']
    body = os.environ['PR_BODY']
    bot = os.environ['BOT_AUTHOR']
except KeyError:
    raise ValueError('PR_NUMBER, PR_BODY and BOT_AUTHOR env vars must be set.')
assert str(int(pr)) == pr
body = body.replace('\\r', '').replace('\\n', '\n')

# Get changelog directory
# TODO: get directory from pyproject.toml
path = os.path.join(os.getenv('CHANGELOG_DIR', 'changelog'))
if not os.path.exists(path):
    raise ValueError(f"CHANGELOG_DIR path '{path}' must exist.")

# Check if bot has control of the changelogs
existing = glob.glob(os.path.join(path, f"{pr}.*.rst"))
for entry in existing:
    author = subprocess.check_output(["git", "log", "-s", "-n1", "--pretty=%an <%ae>", entry], text=True).strip()
    if author != bot:
        print(f"{bot!r} is not modifying any changelogs because")
        print(f"{entry!r} was last modified by {author!r}.")
        exit(0)

# Remove existing
for entry in existing:
    os.remove(entry)
    print(f"Deleted existing changelog: {entry}")

# Search for changelog entries
#     markdown heading of "Changelog [category]"
#     followed by a message in the next ``` code block.
body = re.sub("(<!--.*?-->)", "", body, flags=re.DOTALL)
regex = r"[#]+[ ]*Changelog[ ]*\[([a-zA-Z]+)\][^`]*```((?s:.)*?)```"
results = re.findall(regex, body)

# Count entries per category
count = {}
for category, _ in results:
    inc(count, category)

# Generates (filename, message) pairs
entries = []  # list of changelogs
n = {}  # counter for each category
for category, text in results:

    if count[category] > 1:  # include filename [.COUNTER]
        inc(n, category)
        counter = f'.{n[category]}'
    else:
        counter = ''

    filename = f"{pr}.{category}{counter}.rst"
    text = text.strip().replace('\r', '')
    if len(text) > 0:  # ignore empty messages
        entries += [(filename, text)]

if len(entries) == 0:
    print("No new changelog entries found.")
    exit(0)

for filename, text in entries:
    full_path = os.path.join(path, filename)
    with open(full_path, 'w') as f:
        f.write(text + '\n')
    print(f"Wrote new changelog: {full_path}")
