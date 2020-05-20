import os

import sphinx.util.docutils as sphinx_docutils
from docutils import statemachine
from docutils.io import FileInput
from docutils.parsers.rst import Directive

__all__ = ['MiniGallery']


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


def setup(app):
    # If a minigallery directive has already been registered, do nothing.
    # It's probably that we are on a version of sphinx-gallery with the directive
    if sphinx_docutils.is_directive_registered("minigallery"):
        print("Not using sunpy.util.sphinx.minigallery as a minigallery "
              "directive has already been registered.")
        return

    app.add_directive('minigallery', MiniGallery)
