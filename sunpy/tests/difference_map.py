import sunpy
import matplotlib.colors as colors

f1 = '/Users/schriste/Dropbox/AIA20101016_191231_0193.fits'
f2 = '/Users/schriste/Dropbox/AIA20101016_191222_0193.fits'

map1 = sunpy.Map(f1)
map2 = sunpy.Map(f2)

dmap = map2 - map1

dmap.show(norm = colors.Normalize(-5,5,True))