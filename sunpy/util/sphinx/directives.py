import os
import sys
from io import StringIO

from docutils import nodes, statemachine
from docutils.io import FileInput
from docutils.parsers.rst import Directive, directives

__all__ = ['MiniGallery', 'Generate']


class Generate(Directive):
    """
    Custom directive to include raw output generated using supplied Python code

    This directive is similar to the ``raw`` directive, except that the raw content is generated
    instead of statically provided.  As with the ``raw`` directive, the required argument specifies
    the format of the output (e.g., ``html`` or ``latex``).

    The optional flag ``html_border`` surrounds HTML output with a black border for visual
    separation.
    """
    required_arguments = 1
    optional_arguments = 0
    option_spec = {'html_border': directives.flag}
    has_content = True

    def run(self):
        # Respect the same disabling options as the ``raw`` directive
        if (not self.state.document.settings.raw_enabled
                or not self.state.document.settings.file_insertion_enabled):
            raise self.warning('"%s" directive disabled.' % self.name)

        attributes = {'format': ' '.join(self.arguments[0].lower().split())}

        source, lineno = self.state_machine.get_source_and_line(self.lineno)

        old_stdout, sys.stdout = sys.stdout, StringIO()
        try:
            # Exceute the Python code and capture its output
            exec('\n'.join(self.content))
            text = sys.stdout.getvalue()

            # Wrap HTML output in a black border if requested
            if attributes['format'] == 'html' and 'html_border' in self.options:
                text = f"<div style='border:1px solid black; padding:3px'>{text}</div>"

            # Append the output in the same fashion as the ``raw`` directive
            raw_node = nodes.raw('', text, **attributes)
            raw_node.source, raw_node.line = source, lineno
            return [raw_node]
        except Exception as e:
            message = f"Unable to execute Python code at {os.path.basename(source)}:{lineno}"
            return [nodes.error(None, nodes.paragraph(text=message)),
                    nodes.paragraph(text=str(e))]
        finally:
            sys.stdout = old_stdout


class MiniGallery(Directive):
    """
    Custom directive to insert a mini-gallery

    The required argument is the qualified name of the object.  The mini-gallery will be the subset
    of gallery examples that make use of that object (from that specific namespace).  There can be
    more than one object named, separated by spaces.
    """
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True

    def run(self):
        # Respect the same disabling options as the ``raw`` directive
        if (not self.state.document.settings.raw_enabled
                or not self.state.document.settings.file_insertion_enabled):
            raise self.warning('"%s" directive disabled.' % self.name)

        obj_list = self.arguments[0].split()

        # Concatenate the backreferences file(s)
        lines = []
        for obj in obj_list:
            path = os.path.join(os.getcwd(), 'generated', 'modules', f'{obj}.examples')
            lines += (FileInput(source_path=path).readlines())[5:]  # slice removes heading

        # Append the end for the gallery
        lines += ['\n',
                  '.. raw:: html\n',
                  '\n',
                  '    <div class="sphx-glr-clear"></div>\n']
        text = ''.join(lines)

        include_lines = statemachine.string2lines(text, convert_whitespace=True)
        self.state_machine.insert_input(include_lines, path)

        return []
