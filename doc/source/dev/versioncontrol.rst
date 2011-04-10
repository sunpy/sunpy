Version Control
---------------

Source-code for SunPy is managed using the Bazaar Distributed Version Control 
System (VCS). Code branches are hosted on Launchpad.net, a free project hosting 
website for Open-Source software.

Collaboration
-------------

When multiple people are working on SunPy at the same time, the methods 
described by the Team Collaboration/Distributed Development article should be
used as defined in the Bazar User Guide.

Each developer should has his or her own personal development branch (e.g. 
"john-dev") where all of the main coding is done. When they have finished making
changes and the code has been tested and verified to be working well, the code 
can be merged back into the main trunk. When multiple developers are working on 
SunPy at the same time, care should be taken to avoid merging problems. Each user
should keep a copy of the trunk on their own machine. Each day, the programmer
can use "bzr pull" on the trunk followed by "bzr merge" on their development 
branch to include any recent changes made by other developers.
