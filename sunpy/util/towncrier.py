"""
This module provides Utility functions for working with the towncrier
changelog.

This file is based heavily on towncrier, please see
licenses/TOWNCRIER.rst
"""
import os
from functools import partial

import pkg_resources
from towncrier import (_get_date, append_to_newsfile, find_fragments, get_project_name,
                       get_version, load_config, render_fragments, split_fragments)

__all__ = ["generate_changelog_for_docs"]


def generate_changelog_for_docs(directory, output_filename=None):
    """
    This is a modified version of the `towncrier._main` function with a few
    things disabled.
    """
    print("Updating Changelog...")
    directory = os.path.abspath(directory)
    _join_dir = partial(os.path.join, directory)
    config = load_config(directory)

    print("Loading template...")
    if config["template"] is None:
        template = pkg_resources.resource_string("towncrier", "templates/template.rst").decode("utf8")
    else:
        with open(config["template"], "rb") as tmpl:
            template = tmpl.read().decode("utf8")

    print("Finding news fragments...")

    definitions = config["types"]

    if config.get("directory"):
        base_directory = _join_dir(config["directory"])
        fragment_directory = None

    fragments, fragment_filenames = find_fragments(base_directory, config["sections"],
                                                   fragment_directory, definitions)

    print("Rendering news fragments...")
    fragments = split_fragments(fragments, definitions)
    rendered = render_fragments(
        # The 0th underline is used for the top line
        template,
        config["issue_format"],
        fragments,
        definitions,
        config["underlines"][1:],
        config["wrap"],
    )

    project_version = get_version(_join_dir(config["package_dir"]), config["package"])

    package = config.get("package")
    if package:
        project_name = get_project_name(
            os.path.abspath(_join_dir(config["package_dir"])), package)
    else:
        # Can't determine a project_name, but maybe it is not needed.
        project_name = ""

    project_date = _get_date()

    top_line = config["title_format"].format(
        name=project_name, version=project_version, project_date=project_date)
    top_line += u"\n" + (config["underlines"][0] * len(top_line)) + u"\n"

    print("Writing to newsfile...")
    start_line = config["start_line"]
    if not output_filename:
        output_filename = _join_dir(config["filename"])
    output_filename = os.path.abspath(output_filename)
    append_to_newsfile(directory, output_filename, start_line, top_line, rendered)

    print("Done!")
