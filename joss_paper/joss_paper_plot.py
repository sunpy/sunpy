import matplotlib
import matplotlib.pyplot as plt
import pylab

from astropy import units as u
from astropy.coordinates import SkyCoord

import sunpy.data.sample
import sunpy.map
import sunpy.timeseries as ts
from sunpy.net import hek
from sunpy.time import parse_time

client = hek.HEKClient()

# GOES data for timeseries
goes = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')

flares_hek = client.search(hek.attrs.Time('2011-06-07 00:00', '2011-06-07 23:59'),
                           hek.attrs.FL, hek.attrs.FRM.Name == 'SWPC')

# AIA data for map
my_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
top_right = SkyCoord(1200 * u.arcsec, 0 * u.arcsec, frame=my_map.coordinate_frame)
bottom_left = SkyCoord(500 * u.arcsec, -700 * u.arcsec, frame=my_map.coordinate_frame)
my_submap = my_map.submap(bottom_left, top_right)


# plot figure
fig = plt.figure(figsize=(13, 6))

ax0 = pylab.axes([0.05, 0.09, 0.42, 0.8])
ax0.plot(goes.data['xrsb'], color='r', label=r'1-8 $\mathrm{\AA}$')
ax0.plot(goes.data['xrsa'], color='b', label=r'0.5-4 $\mathrm{\AA}$')

ax0.set_yscale('log')
ax0.set_ylim(1e-9, 1e-3)
ax0.set_ylabel('Watts m$^{-2}$')
ax0.set_xlabel('Time (UT) 2011-06-07')
ax0.set_title('GOES X-ray flux')

ax0.axvline(parse_time(flares_hek[0].get('event_peaktime')
                       ).to_datetime(), ls='dashed', color='grey', label='Flare peak')

ax0.yaxis.grid(True, 'major')
ax0.xaxis.grid(False, 'major')
ax0.legend()

ax0.xaxis.set_tick_params(rotation=30)
ax0.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))
ax0.set_xlim('2011-06-07 00:00', '2011-06-07 23:59')

ax_00 = ax0.twinx()
ax_00.set_yscale("log")
ax_00.set_ylim(1e-9, 1e-3)
ax_00.set_yticks((1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3))
ax_00.set_yticklabels((' ', 'A', 'B', 'C', 'M', 'X'))


ax1 = pylab.axes([0.52, 0.08, 0.5, 0.82], projection=my_submap)
my_submap.plot(clip_interval=(1, 99.95)*u.percent, title='', axes=ax1)
my_submap.draw_grid(axes=ax1)
my_submap.draw_limb(axes=ax1)
cbar_ax = plt.gcf().add_axes([0.97, 0.08, 0.01, 0.82])
plt.colorbar(cax=cbar_ax)


plt.savefig('joss_paper_plot.pdf', dpi=200)
plt.close()
