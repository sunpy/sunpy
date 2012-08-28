from __future__ import absolute_import

import os
import glob

import sunpy
import urllib2

from sunpy.net.util import get_system_filename
from sunpy.util.util import buffered_write, replacement_filename

from sunpy.util.cond_dispatch import ConditionalDispatch, run_cls

class Parent(object):
    _create = ConditionalDispatch()

    @classmethod
    def read(self, filename):
        raise NotImplementedError

    @classmethod
    def from_glob(cls, pattern):
        """ Read out files using glob (e.g., ~/BIR_2011*) pattern. Returns 
        list of CallistoSpectrograms made from all matched files.
        """        
        return map(cls.read, glob.glob(pattern))
    
    @classmethod
    def from_single_glob(cls, singlepattern):
        """ Read out a single file using glob (e.g., ~/BIR_2011*) pattern.
        If more than one file matches the pattern, raise ValueError.
        """
        matches = glob.glob(os.path.expanduser(singlepattern))
        if len(matches) != 1:
            raise ValueError("Invalid number of matches: %d" % len(matches))
        return cls.read(matches[0])
    
    @classmethod
    def from_files(cls, filenames):
        """ Return list of CallistoSpectrogram read from given list of
        filenames. """
        filenames = map(os.path.expanduser, filenames)
        return map(cls.read, filenames)
    
    @classmethod
    def from_file(cls, filename):
        """ Return CallistoSpectrogram from FITS file. """
        filename = os.path.expanduser(filename)
        return cls.read(filename)
    
    @classmethod
    def from_dir(cls, directory):
        """ Return list that contains all files in the directory read in. """
        directory = os.path.expanduser(directory)
        return map(cls.read,
            (os.path.join(directory, elem) for elem in os.listdir(directory))
        )
    
    @classmethod
    def from_url(cls, url):
        """ Return CallistoSpectrogram read from URL.
        
        Parameters
        ----------
        url : str
            URL to retrieve the data from
        """
        default_dir = sunpy.config.get("downloads", "download_dir")
        opn = urllib2.urlopen(url)
        try:
            name = get_system_filename(opn, url)
            path = os.path.join(default_dir, name)
            if os.path.exists(path):
                path = replacement_filename(path)

            with open(path, 'wb') as fd:
                buffered_write(opn, fd, 9096)

            return cls.read(path)
        finally:
            opn.close()


Parent._create.add(
    run_cls('from_file'),
    lambda cls, filename: os.path.isfile(os.path.expanduser(filename)),
    [type, basestring], check=False
)
Parent._create.add(
# pylint: disable=W0108
# The lambda is necessary because introspection is peformed on the
# argspec of the function.
    run_cls('from_dir'),
    lambda cls, directory: os.path.isdir(os.path.expanduser(directory)),
    [type, basestring], check=False
)
# If it is not a kwarg and only one matches, do not return a list.
Parent._create.add(
    run_cls('from_single_glob'),
    lambda cls, singlepattern: ('*' in singlepattern and
                           len(glob.glob(
                               os.path.expanduser(singlepattern))) == 1),
    [type, basestring], check=False
)
# This case only gets executed under the condition that the previous one wasn't.
# This is either because more than one file matched, or because the user
# explicitely used pattern=, in both cases we want a list.
Parent._create.add(
    run_cls('from_glob'),
    lambda cls, pattern: '*' in pattern and glob.glob(
        os.path.expanduser(pattern)
        ),
    [type, basestring], check=False
)
Parent._create.add(
    run_cls('from_files'),
    lambda cls, filenames: True,
    types=[type, list], check=False
)
Parent._create.add(
    run_cls('from_url'),
    lambda cls, url: True,
    types=[type, basestring], check=False
)
