# -*- coding: utf-8 -*-
# Author: David Perez-Suarez <dps.helio-?-gmail.com>

from __future__ import absolute_import

import os
import datetime

import numpy as np

from sunpy.util.cond_dispatch import ConditionalDispatch
from sunpy.spectra.spectrogram import LinearTimeSpectrogram, REFERENCE, get_day

__all__ = ['SWavesSpectrogram']

#def query(start, end, spacecraft=None, band=None):
#    start = datetime.datetime(start.year, start.month, start.day)
#    end = datetime.datetime(end.year, end.month, end.day)
#    if spacecraft is None:
#        spacecraft = ['STEREO_A','STEREO_B']
#    elif spacecraft.upper() != 'STEREO_B':
#        spacecraft = ['STEREO_A']
#    else:
#        spacecraft = ['STEREO_B']
#    if band is None:
#        band = 'All'
#        min_wave = 0.00001
#        max_wave = 0.016
#    elif band.lower() != 'low':
#        band = 'high' # range 0.125 - 16 MHz
#        min_wave = 0.00017
#        max_wave = 0.016
#    elif band.lower() == 'low': # range 10 - 160 KHz
#        min_wave = 0.00001
#        max_wave = 0.00012
#    print "Getting data from %s satellite, band %s" % (spacecraft,band)
#    client = vso.InteractiveVSOClient()
#    result = client.query_legacy(tstart=start,\
#                                     tend=end,instrument='SWAVES',source=spacecraft.upper(),\
#                                     min_wave=min_wave,max_wave=max_wave,unit_wave='GHz')
#    if len(result) > 0:
#        listfiles = client.get(result).wait()


class SWavesSpectrogram(LinearTimeSpectrogram):
    _create = ConditionalDispatch.from_existing(LinearTimeSpectrogram._create)
    create = classmethod(_create.wrapper())
    COPY_PROPERTIES = LinearTimeSpectrogram.COPY_PROPERTIES + [
        ('bg', REFERENCE)
    ]

    @staticmethod
    def swavesfile_to_date(filename):
        _, name = os.path.split(filename)
        date = name.split('_')[2]
        return datetime.datetime(
            int(date[0:4]), int(date[4:6]), int(date[6:])
        )

    @classmethod
    def read(cls, filename, **kwargs):
        """Read in FITS file and return a new SWavesSpectrogram. """
        data = np.genfromtxt(filename, skip_header=2)
        time_axis = data[:, 0] * 60.
        data = data[:, 1:].transpose()
        header = np.genfromtxt(filename, skip_footer=time_axis.size)
        freq_axis = header[0, :]
        bg = header[1, :]
        start = cls.swavesfile_to_date(filename)
        end = start + datetime.timedelta(seconds=time_axis[-1])
        t_delt = 60.
        t_init = (start - get_day(start)).seconds
        content = ''
        t_label = 'Time [UT]'
        f_label = 'Frequency [KHz]'

        freq_axis = freq_axis[::-1]
        data = data[::-1, :]

        return cls(data, time_axis, freq_axis, start, end, t_init, t_delt,
            t_label, f_label, content, bg)


    def __init__(self, data, time_axis, freq_axis, start, end,
            t_init, t_delt, t_label, f_label, content, bg):
        # Because of how object creation works, there is no avoiding
        # unused arguments in this case.
        # pylint: disable=W0613

        super(SWavesSpectrogram, self).__init__(
            data, time_axis, freq_axis, start, end,
            t_init, t_delt, t_label, f_label,
            content, set(["SWAVES"])
        )
        self.bg = bg


SWavesSpectrogram.create.im_func.__doc__ = (
    """ Create SWavesSpectrogram from given input dispatching to the
    appropriate from_* function.

Possible signatures:

""" + SWavesSpectrogram._create.generate_docs()
)

if __name__ == "__main__":
    opn = SWavesSpectrogram.read("/home/florian/swaves_average_20120705_a_hfr.dat")
    opn.plot(min_=0, linear=False).show()
    print "Press return to exit"
    raw_input()
