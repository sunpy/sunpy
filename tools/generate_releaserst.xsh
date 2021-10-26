#!/usr/bin/env xonsh
"""
Re-write release.rst for SunPy


Usage:
    generate_releaserst.xsh <prev_version> [<prev_tag>] [--project-name=<project-name>] [--author-sort=<author-sort>] [--show-commit-count] [--pretty-project-name=<pretty-project-name>] [--repo=<repo>] [--pat=<token>]

Options:
    prev_version                                   The PyPI release name of the previous release (should not start with v)
    prev_tag                                       The tag name for the previous release, if not specified will be v<prev_version>
    --project-name=<project-name>                  The project name on PyPI [default: sunpy]
    --author-sort=<author-sort>                    How to sort the contributors in the list. Can be alphabet or numeric. [default: alphabet]
    --show-commit-count                            Show number of commits next to the contributors [default: false]
    --pretty-project-name=<pretty-project-name>    The project name to use in the printed output [default: <project-name>]
    --repo=<repo>                                  The GitHub repository name, will default to <project-name>/<project-name>
    --pat=<token>                                  A GitHub Personal Access Token
"""
# The GitHub stuff is lovingly stolen from astropy-procedures

import os
import json
import netrc
import getpass
import argparse
import datetime
import warnings
from functools import partial

import docopt
import requests

args = docopt.docopt(__doc__, argv=$ARGS[1:], version="sunpy")

ISO_FORMAT = "%Y-%m-%dT%H:%M:%SZ"


def do_graphql_query(owner, repository, etype, token):
    if etype == "pull":
        query_template = """
{{
  repository(owner: "{owner}", name: "{repository}") {{
    pullRequests(first:100, orderBy: {{direction: ASC, field: CREATED_AT}}, baseRefName: "main", states: MERGED{after}) {{
      edges {{
        node {{
          title
          number
          createdAt
          updatedAt
          mergedAt
        }}
        cursor
      }}
    }}
  }}
}}"""
    elif etype == "issue":
        query_template = """
{{
  repository(owner: "{owner}", name: "{repository}") {{
    issues(first:100, orderBy: {{direction: ASC, field: CREATED_AT}}, states: CLOSED{after}) {{
      edges {{
        node {{
          title
          number
          createdAt
          updatedAt
          closedAt
        }}
        cursor
      }}
    }}
  }}
}}"""
    else:
        raise ValueError()


    results = {}
    cursor = ''
    headers = {"Authorization": f"Bearer {token}"}
    while True:

        if not cursor:
            after = ''
        else:
            after = f', after:"{cursor}"'

        query = query_template.format(owner=owner, repository=repository, after=after)
        request = requests.post('https://api.github.com/graphql', json={'query': query}, headers=headers)

        if request.status_code != 200:
            raise Exception(f"Query failed {request.status_code}")

        entries = request.json()
        if 'errors' in entries:
            print(entries)
        if etype == "pull":
            entries = entries['data']['repository']['pullRequests']['edges']
        if etype == "issue":
            entries = entries['data']['repository']['issues']['edges']

        for entry in entries:
            item = entry['node']
            cursor = entry['cursor']
            res = {}

            # Convert times to datetime objects
            for key, value in item.items():
                if key.endswith('At'):
                    value = datetime.datetime.strptime(value, ISO_FORMAT)
                res[key] = value

            results[item['number']] = res

        if len(entries) < 100:
            break

    return results

def filter_between_dates(since, upto, key, item):
    number, info = item
    time = info[key]
    return since < time and time < upto

def count_issues_since(since, upto, repo, token, verbose):
    owner, repository = repo.split('/')
    return len(list(
        filter(partial(filter_between_dates, since, upto, "closedAt"),
               do_graphql_query(owner, repository, "issue", token).items())
        ))

def count_prs_since(since, upto, repo, token, verbose):
    owner, repository = repo.split('/')
    return len(list(
        filter(partial(filter_between_dates, since, upto, "mergedAt"),
               do_graphql_query(owner, repository, "pull", token).items())
    ))


def get_datetime_of_pypi_version(pkg, version):
    resp = requests.get(f"https://pypi.org/pypi/{pkg}/json")
    j = resp.json()

    datestr = j['releases'][version][0]['upload_time']

    return datetime.datetime.strptime(datestr, "%Y-%m-%dT%H:%M:%S")



# Handle input

prev_version = args['<prev_version>']
if prev_version.startswith("v"):
    raise ValueError("prev_version should not start with v")
prev_tag = args['<prev_tag>'] if args['<prev_tag>'] else "v"+args['<prev_version>']

commit_count = args['--show-commit-count']

author_sort = args['--author-sort']
if author_sort not in ("alphabet", "numeric"):
    raise ValueError("--author_sort should be one of 'alphabet' or 'numeric'")


# Parse Git trickery.
flags = "-s"
# Sort numerically if the option is specified
flags = flags + "n" if author_sort == "numeric" else flags
current_log = $(git shortlog @(flags) --no-merges @(prev_tag+"..HEAD"))

# Get all the authors for all releases up to the previous release
prev = {a.split('\t')[1] for a in $(git shortlog -ns --no-merges @($(git rev-list --max-parents=0 HEAD).split("\n")[0]+".."+prev_tag)).split('\n')[:-1]}

# Get all authors from the previous release to this one
current = {a.split('\t')[1] for a in current_log.split('\n')[:-1]}

new = current.difference(prev)

gitcount = $(git rev-list HEAD @("^"+prev_tag) --count)
ncommits = int(gitcount.strip())
npeople = len(current)
nnew = len(new)

# Reformat the log

lines = current_log.split('\n')[:-1]

shortlog = []
for i, line in enumerate(lines):
    if commit_count:
        outl = line
    else:
        outl = line.split('\t')[1]
    if any([a in line for a in new]):
        outl += '  *'
    shortlog.append(outl)

shortlog = list(map(lambda x: '-  ' + x, shortlog))

# Get PR info

since = get_datetime_of_pypi_version(args['--project-name'], prev_version)
upto = datetime.datetime.fromisoformat($(git show -s --format=%cI HEAD).strip()).astimezone().replace(tzinfo=None)

verbose = False
repo = f"{args['--project-name']}/{args['--project-name']}" if not args['--repo'] else args['--repo']
icnt = count_issues_since(since, upto, repo, args['--pat'], verbose=verbose)
prcnt = count_prs_since(since, upto, repo, args['--pat'], verbose=verbose)

# Build output
output = '\n'.join(shortlog)


pretty_project_name = args["--pretty-project-name"] if args["--pretty-project-name"] else args["--project-name"]

print()
print(f"This release of {pretty_project_name} contains {ncommits} commits in {prcnt} merged pull requests closing {icnt} issues from {npeople} people, {nnew} of which are first-time contributors to {pretty_project_name}.")
print()
print(f"* {ncommits} commits have been added since {prev_version[:3]}")
print(f"* {icnt} issues have been closed since {prev_version[:3]}")
print(f"* {prcnt} pull requests have been merged since {prev_version[:3]}")
print(f"* {npeople} people have contributed since {prev_version[:3]}")
print(f"* {nnew} of which are new contributors")
print()
print("The people who have contributed to the code for this release are:")
print()
print(output)
print()
print(f"Where a * indicates that this release contains their first contribution to {pretty_project_name}.")
