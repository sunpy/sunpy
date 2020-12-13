import tempfile

from docutils import statemachine
from docutils.parsers.rst import Directive

from sunpy.util.towncrier import generate_changelog_for_docs

__all__ = ['ChangeLog']


class ChangeLog(Directive):
    """
    Render the changelog for the current commit using towncrier.

    This directive renders all the towncrier newsfiles into your current
    documentation, this can be used to keep a rendered version of the changelog
    since your last release in your documentation.

    The directive takes one argument which is the location of your
    ``pyproject.toml`` file (towncrier configuration) relative to the
    ``conf.py`` file *not* the file in which the directive is located.
    If this argument is not specified it defaults to `"../"`.

    Examples
    --------

    .. code-block:: rst

        .. changelog::
    """
    required_arguments = 0
    optional_arguments = 1
    final_argument_whitespace = True

    def run(self):
        config_path = self.arguments[0] or "../"
        output_file = tempfile.mkstemp()[1]
        generate_changelog_for_docs(config_path, output_filename=output_file)
        with open(output_file) as fobj:
            include_lines = statemachine.string2lines(fobj.read(), convert_whitespace=True)

        self.state_machine.insert_input(include_lines, "")

        return []


class DummyChangelog(ChangeLog):
    def run(self):
        return []


def setup(app):
    app.add_directive('changelog', ChangeLog)

    return {'parallel_read_safe': True, 'parallel_write_safe': True}
