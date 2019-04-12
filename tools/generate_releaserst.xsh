#!/usr/bin/env xonsh
"""
Re-write release.rst for SunPy


Usage:
    generate_releasemd.xsh <prev_version> [<prev_tag>] [--commit-count] [--project-name=<project-name>] [--pretty-project-name=<pretty-project-name>] [--repo=<repo>] [--auth]

Options:
    prev_version                                   The PyPI release name of the previous release
    prev_tag                                       The tag name for the previous release, if not specified will be v<prev_version>
    --commit-count                                 Generate commit count stats
    --project-name=<project-name>                  The project name on PyPI [default: sunpy]
    --pretty-project-name=<pretty-project-name>    The project name to use in the printed output [default: SunPy]
    --repo=<repo>                                  The GitHub repository name, will default to <project-name>/<project-name>
    --auth                                         Prompt for GH auth
"""
# The GitHub stuff is lovingly stolen from astropy-procedures

import netrc
import getpass
import warnings
import argparse
import os
import json
import datetime

import docopt
import requests

args = docopt.docopt(__doc__, argv=$ARGS[1:], version="sunpy")

GH_API_BASE_URL = 'https://api.github.com'
ISO_FORMAT = "%Y-%m-%dT%H:%M:%SZ"

def get_credentials(username=None, password=None):

    try:
        my_netrc = netrc.netrc()
    except Exception as e:
        print(e)
    else:
        auth = my_netrc.authenticators("api.github.com")
        if auth:
            response = 'NONE'  # to allow enter to be default Y
            while response.lower() not in ('y', 'n', ''):
                warnings.warn('Using the following GitHub credentials from '
                              '~/.netrc: {0}/{1}'.format(auth[0], '*' * 8))
                response = input(
                    'Use these credentials (if not you will be prompted '
                    'for new credentials)? [Y/n] ')
            if response.lower() == 'y' or response == '':
                username = auth[0]
                password = auth[2]

    if not (username or password):
        warnings.warn("Enter your GitHub username and password so that API "
                      "requests aren't as severely rate-limited...")
        username = input('Username: ')
        password = getpass.getpass('Password: ')
    elif not password:
        warnings.warn("Enter your GitHub password so that API "
                      "requests aren't as severely rate-limited...")
        password = getpass.getpass('Password: ')

    return username, password

auth = get_credentials() if args['--auth'] else None


def paginate_list_request(req, verbose=False, auth=None):
    elems = []
    currreq = req
    i = 1

    while 'next' in currreq.links:
        elems.extend(currreq.json())

        i += 1
        if verbose:
            print('Doing request', i, 'of', currreq.links['last']['url'].split('page=')[-1])
        currreq = requests.get(currreq.links['next']['url'], auth=auth)

    elems.extend(currreq.json())
    return elems


def count_issues_since(dt, repo, auth=None, verbose=True, cacheto=None):
    if cacheto and os.path.exists(cacheto):
        with open(cacheto) as f:
            isslst = json.load(f)
    else:
        url = GH_API_BASE_URL + '/repos/' + repo + '/issues?per_page=100&state=all'

        req = requests.get(url, auth=auth)
        if not req.ok:
            msg = 'Failed to access github API for repo using url {}. {}: {}: {}'
            raise requests.HTTPError(msg.format(url, req.status_code, req.reason, req.text))
        isslst = paginate_list_request(req, verbose, auth=auth)
        if cacheto:
            with open(cacheto, 'w') as f:
                json.dump(isslst, f)

    nopened = nclosed = 0

    for entry in isslst:
        if not isinstance(entry, dict):
            continue
        createddt = datetime.datetime.strptime(entry['created_at'],  ISO_FORMAT)
        if createddt > dt:
            nopened += 1

        if entry['closed_at']:
            closeddt = datetime.datetime.strptime(entry['closed_at'],  ISO_FORMAT)
            if closeddt > dt:
                nclosed += 1

    return {'opened': nopened, 'closed': nclosed}


def count_prs_since(dt, repo, auth=None, verbose=True, cacheto=None):
    if cacheto and os.path.exists(cacheto):
        with open(cacheto) as f:
            prlst = json.load(f)
    else:
        url = GH_API_BASE_URL + '/repos/' + repo + '/pulls?per_page=100&state=all'

        req = requests.get(url, auth=auth)
        prlst = paginate_list_request(req, verbose, auth=auth)
        if cacheto:
            with open(cacheto, 'w') as f:
                json.dump(prlst, f)


    nopened = nclosed = 0

    usersopened = []
    usersclosed = []

    for entry in prlst:
        if not isinstance(entry, dict):
            continue
        createddt = datetime.datetime.strptime(entry['created_at'],  ISO_FORMAT)
        if createddt > dt:
            nopened += 1
            user = entry['user']
            if user is not None:
                usersopened.append(user['id'])

        if entry['merged_at']:
            closeddt = datetime.datetime.strptime(entry['merged_at'],  ISO_FORMAT)
            if closeddt > dt:
                nclosed += 1
                user = entry['user']
                if user is not None:
                    usersclosed.append(user['id'])

    return {'opened': nopened, 'merged': nclosed,
            'usersopened': len(set(usersopened)),
            'usersmerged': len(set(usersclosed))}


def get_datetime_of_pypi_version(pkg, version):
    resp = requests.get(f"https://pypi.org/pypi/{pkg}/json")
    j = resp.json()

    datestr = j['releases'][version][0]['upload_time']

    return datetime.datetime.strptime(datestr, "%Y-%m-%dT%H:%M:%S")



# Handle input

prev_version = args['<prev_version>']
prev_tag = args['<prev_tag>'] if args['<prev_tag>'] else "v"+args['<prev_version>']

commit_count = args['--commit-count']


# Parse Git trickery.

current_log = $(git shortlog -ns --no-merges @(prev_tag+"..HEAD"))

# Get all the authors for all releases up to the previous release
prev = {a.split('\t')[1] for a in $(git shortlog -ns --no-merges @($(git rev-list --max-parents=0 HEAD).split("\n")[0]+".."+prev_tag)).split('\n')[:-1]}

# Get all authors from the previous release to this one
current = {a.split('\t')[1] for a in current_log.split('\n')[:-1]}

new = current.difference(prev)

ncommits = int($(git rev-list HEAD @("^"+prev_tag) --count).strip())
npeople = len(current)
nnew = len(new)

# Reformat the log

lines = current_log.split('\n')[:-1]

shortlog = []
for i, line in enumerate(lines):
    if commit_count:
        outl = line
    else:
        outl = '    ' + line.split('\t')[1]
    if any([a in line for a in new]):
        outl += '  *'
    shortlog.append(outl)


# Get PR info

pkgdt = get_datetime_of_pypi_version(args['--project-name'], prev_version)

icache = 'issues.json'
prcache = 'prs.json'

verbose = False
repo = f"{args['--project-name']}/{args['--project-name']}" if not args['--repo'] else args['--repo']
icnt = count_issues_since(pkgdt, repo, auth=auth, verbose=verbose, cacheto=icache)
prcnt = count_prs_since(pkgdt, repo, auth=auth, verbose=verbose, cacheto=prcache)

# Build output
output = '\n'.join(shortlog)


pretty_project_name = args["--pretty-project-name"] if args["--pretty-project-name"] else args["--project-name"]
print(f"This release of {pretty_project_name} contains {ncommits} commits in {prcnt['merged']} merged pull requests closing {icnt['closed']} issues from {npeople} people, {nnew} of which are first time contributors to {pretty_project_name}.")
print()
print("The people who have contributed to the code for this release are:")
print()
print(output)
print()
print(f"Where a * indicates their first contribution to {pretty_project_name}.")
