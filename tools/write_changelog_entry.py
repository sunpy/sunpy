import glob
import os
import re


def inc(obj, key):
    """Increment a dictionary[key] value."""
    try:
        obj[key] += 1
    except KeyError:
        obj[key] = 1


# Load PR number and comment
try:
    pr = os.environ['PR_NUMBER']
    body = os.environ['PR_BODY']
except KeyError:
    raise ValueError('PR_NUMBER and PR_BODY env vars must be set.')
assert str(int(pr)) == pr
body = body.replace('\\r', '').replace('\\n', '\n')

# Get changelog directory
# TODO: get directory from pyproject.toml
path = os.path.join(os.getenv('CHANGELOG_DIR', 'changelog'))
if not os.path.exists(path):
    raise ValueError(f"CHANGELOG_DIR path '{path}' must exist.")

# Search for changelog entries
#     markdown heading of "Changelog [category]"
#     followed by a message in the next ``` code block.
body = re.sub("(<!--.*?-->)", "", body, flags=re.DOTALL)
regex = r"[#]+[ ]*Changelog[ ]*\[([a-zA-Z]+)\]*[^`]*```([^`]*)```"
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
    text = text.strip()
    if len(text) > 0:  # ignore empty messages
        entries += [(filename, text)]

for existing in glob.glob(os.path.join(path, f"{pr}.*.rst")):
    os.remove(existing)
    print(f"Deleted exiting changelog: {existing}")

if len(entries) == 0:
    raise ValueError("No new changelog entries found.")

for filename, text in entries:
    full_path = os.path.join(path, filename)
    with open(full_path, 'w') as f:
        f.write(text + '\n')
    print(f"Wrote new changelog: {full_path}")
