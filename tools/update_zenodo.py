import json
import subprocess


def remove_initials(name):
    # Remove initials for a string name
    # Assumes names/initials are all separated by a single space
    new_name = []
    for n in name.split(' '):
        if len(n) == 2 and n[1] == '.':
            continue
        new_name.append(n)
    return ' '.join(new_name)


authors = subprocess.check_output(['git shortlog -s -n'], shell=True)
authors = authors.decode('utf-8')
# 1 author per line
authors = authors.split('\n')[:-1]
# Use tab between number of commits and name to get name
authors = [author.split('\t')[1] for author in authors]
# Remove initials
authors = [remove_initials(auth) for auth in authors]
# List of authors to ignore because they are bots
manual_ignore = ['codetriage-readme-bot']
authors = [auth for auth in authors if auth not in manual_ignore]

# Get list of current authors in zenodo.yml
with open('.zenodo.json') as zenodo_file:
    data = json.load(zenodo_file)

creators = data['creators']
already_auth = [auth['name'] for auth in creators]
# Remove any initials
already_auth = [remove_initials(auth) for auth in already_auth]

new_creators = []
# Loop through all the current authors
print('New contributors since the last update of zenodo.json:')
for author in authors:
    # If already in .zenodo.json, take the entry to preserve ORCID and affiliation
    if author in already_auth:
        new_creators.append(creators[already_auth.index(author)])
    else:
        print(author)
        new_creators.append({'name': author})

# Add in anyone currently in .zenodo.json, but who hasn't committed to the repository
print('\nContributors with no code commits:')
non_commitors = list(set(already_auth) - set(authors))
non_commitors.sort()
for author in non_commitors:
    print(author)
    new_creators.append(creators[already_auth.index(author)])

data['creators'] = new_creators
with open('.zenodo.json', 'w') as zenodo_file:
    json.dump(data, zenodo_file, indent='    ', ensure_ascii=False)
