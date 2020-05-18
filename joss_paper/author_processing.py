import io
import copy
import json
from pathlib import Path

from ruamel.yaml import YAML

yaml = YAML()

with open("metadata.yml") as fd:
    data = yaml.load(fd)

authors = data['authors']

# Convert to a zenodo JSON

creators = []

for author in authors:
    creators.append(copy.deepcopy(author))
    if isinstance(author.get('affiliation', ''), list):
        new_affiliation = ''
        for i, affil in enumerate((author['affiliation'])):
            new_affiliation += f"{i + 1} - {affil} "
        creators[-1]['affiliation'] = new_affiliation

with open(Path("../.zenodo.json").resolve(), "r") as fobj:
    zenodo = json.load(fobj)

zenodo['creators'] = creators

with open(Path("../.zenodo.json").resolve(), "w") as fobj:
    json.dump(zenodo, fobj, indent=4, ensure_ascii=False)


# Convert to JOSS metadata

# Create a list of all affiliations and make the affiliation key a list (and always present)
affiliations = []
for auth in authors:
    affil = auth.get('affiliation', "None")
    if affil is not None:
        if not isinstance(affil, list):
            auth['affiliation'] = [affil]

        affiliations += auth['affiliation']

affiliations = set(affiliations)

affiliation_index = {}
for auth in authors:
    new_affils = []
    if auth.get('affiliation') is not None:
        for affil in auth['affiliation']:
            index = affiliation_index.get(affil, len(affiliation_index) + 1)
            affiliation_index[affil] = index
            new_affils.append(index)
    auth['affiliation'] = ', '.join(map(str, new_affils))


# Build the affiliations section from the index
affiliations = []
for name, index in affiliation_index.items():
    affiliations.append({'name': name, 'index': index})

data['affiliations'] = affiliations


# Dump the new yaml to a string
s = io.StringIO()
yaml.default_flow_style = False
yaml.width = 1000
yaml.allow_unicode = True
yaml.dump(data, s)
s.seek(0)


# Overwrite the paper markdown with the new metadata
with open("paper.md") as fd:
    lines = fd.readlines()

start_line = lines.index("---\n") + 1
end_line = lines.index("...\n")

out_lines = lines[:start_line]
out_lines += s.readlines()
out_lines += lines[end_line:]

with open("paper.md", "w") as fd:
    fd.writelines(out_lines)
