import io
# import yaml
from ruamel.yaml import YAML
yaml = YAML()

with open("paper.md") as fd:
    lines = fd.readlines()
    start_line = lines.index("---\n")+1
    end_line = lines.index("...\n")

with open("metadata.yml") as fd:
    data = yaml.load(fd)

authors = data['authors']

# Create a list of all affiliations and make the affiliation key a list (and always present)
affiliations = []
for auth in authors:
    affil = auth.get('affiliation', "The SunPy Developers")
    if affil is not None:
        if not isinstance(affil, list):
            auth['affiliation'] = [affil]

        affiliations += auth['affiliation']

affiliations = set(affiliations)

affiliation_indexes = {"The SunPy Developers": 1}
for auth in authors:
    new_affils = []
    if auth.get('affiliation') is not None:
        for affil in auth['affiliation']:
            index = affiliation_indexes.get(affil, len(affiliation_indexes)+1)
            affiliation_indexes[affil] = index
            new_affils.append(index)
        auth['affiliation'] = ', '.join(map(str, new_affils))


affiliations = []
for name, index in affiliation_indexes.items():
    affiliations.append({'name': name, 'index': index})

data['affiliations'] = affiliations


s = io.StringIO()
yaml.default_flow_style = False
yaml.width = 1000
yaml.allow_unicode = True
yaml.dump(data, s)
s.seek(0)

out_lines = lines[:start_line]
out_lines += s.readlines()
out_lines += lines[end_line:]

with open("paper.md", "w") as fd:
    fd.writelines(out_lines)
