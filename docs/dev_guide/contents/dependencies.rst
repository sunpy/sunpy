.. _grouping_dependencies:

*************************
Extra Packages and Groups
*************************

sunpy has a number of optional dependencies that are not required to run the core functionality of the package.
These dependencies are grouped into "extras" that can be installed with the package.

We have three categories of extras:

* subpackages : These extra groups are named after each subpackage and list the requirements that enable the core functionality of that subpackage.
* optional : These extra groups are named after a package or feature they enable. These are not needed to import sunpy or any of its subpackages, but are needed to use some of the functionality.
  Subpackages will emit a clear import error message if the needed package is not installed.
* groupings : These are the groups that are used to install multiple optional dependencies at once.
  These range from "all" which includes every optional dependency to "core" which only includes the dependencies that are needed to run the core functionality of sunpy.
  This also includes the "test", "docs" and "dev" groups.
