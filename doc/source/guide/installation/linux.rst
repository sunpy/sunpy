=====
Linux
=====
Installation instructions for Linux.

Ubuntu
^^^^^^
To begin, install the pre-requisites for SunPy using :command:`apt-get`: ::

    sudo apt-get install python-numpy python-matplotlib python-pyfits python-qt4 python-scipy python-suds git-core ipython 

The ``ipython`` package in the above list installs the `IPython enhanced console 
<http://ipython.scipy.org/moin/>`_ and is optional but recommended.

Next, use Git to download a copy of the latest version of SunPy: ::

    git clone git://git@github.com/sunpy/sunpy.git

Done! To see if everything went okay, start a Python session and try importing
SunPy:

>>> import sunpy
>>> sunpy.Map(sunpy.AIA_171_IMAGE).show()

Testing your installation
^^^^^^^^^^^^^^^^^^^^^^^^^

Now you can test your installation. Open a new Python shell by typing 
:command:`python` in ``Command Prompt``, and type these commands: ::

>>> import sunpy
>>> sunpy.Map(sunpy.AIA_171_IMAGE).show()