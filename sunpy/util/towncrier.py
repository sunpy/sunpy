"""
This module provides Utility functions for working with the towncrier
changelog.

This file is based heavily on towncrier, please see
licenses/TOWNCRIER.rst
"""
import os
from functools import partial

import pkg_resources
from towncrier._builder import find_fragments, render_fragments, split_fragments
from towncrier._project import get_project_name, get_version
from towncrier._settings import load_config
from towncrier._writer import append_to_newsfile
from towncrier.build import _get_date

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
    if not config:
        raise FileNotFoundError(f"Could not locate the towncrier config file at path {directory}.")

    print("Loading template...")
    if config["template"] is None:
        template = pkg_resources.resource_string(
            "towncrier", "templates/template.rst").decode("utf8")
    else:
        with open(config["template"], "rb") as tmpl:
            template = tmpl.read().decode("utf8")

    print("Finding news fragments...")

    definitions = config["types"]

    if config.get("directory"):
        base_directory = _join_dir(config["directory"])
        fragment_directory = None

    fragments, _ = find_fragments(base_directory, config["sections"],
                                  fragment_directory, definitions)

    print("Rendering news fragments...")
    fragments = split_fragments(fragments, definitions)
    project_version = config.get('version')
    project_version = get_version(
        os.path.join(base_directory, config["package_dir"]), config["package"]
    ).strip()
    project_name = config.get('name')
    if not project_name:
        package = config.get("package")
        if package:
            project_name = get_project_name(
                os.path.abspath(os.path.join(base_directory, config["package_dir"])),
                package,
            )
        else:
            # Can't determine a project_name, but maybe it is not needed.
            project_name = ""
    project_date = _get_date().strip()
    if config["title_format"]:
        top_line = config["title_format"].format(
            name=project_name, version=project_version, project_date=project_date
        )
    else:
        top_line = ""

    rendered = render_fragments(
        # The 0th underline is used for the top line
        template,
        config["issue_format"],
        top_line,
        fragments,
        definitions,
        config["underlines"][1:],
        config["wrap"],
        {"name": project_name, "version": project_version, "date": project_date},
        top_underline=config["underlines"][0],
        all_bullets=config["all_bullets"],
    )

    print("Writing to newsfile...")
    start_line = config["start_string"]
    if not output_filename:
        output_filename = _join_dir(config["filename"])
    output_filename = os.path.abspath(output_filename)
    append_to_newsfile(directory, output_filename, start_line, top_line, rendered)
    print("Done!")
