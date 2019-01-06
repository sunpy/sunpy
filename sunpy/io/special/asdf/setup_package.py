import os
from pathlib import Path


def get_package_data():
    # Installs the schema files
    schemas = []
    root = str(Path.home().joinpath(Path(__file__).parent, 'schemas'))
    for node, dirs, files in os.walk(root):
        for fname in files:
            if fname.endswith('.yaml'):
                schemas.append(
                    os.path.relpath(str(Path.home().joinpath(node, fname)), root))

    # In the package directory, install to the subdirectory 'schemas'
    schemas = [str(Path.home().joinpath('schemas', s)) for s in schemas]

    return {'sunpy.io.special.asdf': schemas}
