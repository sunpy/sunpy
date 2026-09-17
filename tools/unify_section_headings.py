#!/usr/bin/env python

"""
Find the sets of characters used in RST section headers, and replace with a standard sequence.
"""

import re
import shutil
import tempfile
from pathlib import Path

HEADER_CHAR_LEVELS = '*=-^"+:~'


def replace_header_chars(filename):
    header_chars = ""
    level = -1

    rx = re.compile(r"^(\S)\1{3,}\s*$")
    with open(filename) as infile, tempfile.NamedTemporaryFile(
        suffix=".rst", delete=False
    ) as outfile:
        for i, line in enumerate(infile):
            match = rx.match(line)
            if match:
                char = match.group(1)
                if char in header_chars:
                    # Existing header char: go back out to that level
                    new_level = header_chars.index(char)
                    if new_level > level + 1:
                        # doc tries to jump up too many levels
                        raise ValueError(
                            f"ERROR misorder new_level={new_level} level={level} "
                            f"char={char} header_chars={header_chars} on line {i}"
                        )
                    level = new_level
                    print(f"s/{char}/{HEADER_CHAR_LEVELS[level]}/")
                else:
                    # New header char - create a deeper level
                    if level == len(header_chars) - 1:
                        header_chars += char
                        level += 1

                        print(f"s/{char}/{HEADER_CHAR_LEVELS[level]}/")
                    else:
                        # we're trying to create a new level,
                        # but we're not at the current deepest level
                        raise ValueError(
                            f"ERROR misorder {char} at level {level} from {header_chars} on line {i}"
                        )
                outfile.write(line.replace(char, HEADER_CHAR_LEVELS[level]).encode())
            else:
                outfile.write(line.encode())

    print(f"{outfile.name} -> {filename}")
    shutil.move(outfile.name, filename)

    return header_chars


if __name__ == "__main__":
    here = Path(__file__).parent.absolute()
    rst_files = here.glob("../docs/**/*rst")
    [replace_header_chars(filename) for filename in rst_files]
