# Tomas Meszaros <exo@tty.sk>

from sunpy.hypermap.sources.iris import Parser

def make_coord_system():
    p = Parser("./iris_sample_data.fits")
    f = p.get_coordinate_system("Coordinate Test System")

    print f.name
    print len(f.name)*"="
    for i in f.frames:
        print "    system: %s" % i.system
        print "  num_axes: %d" % i.num_axes
        print "axes_names: %s" % i.axes_names
        print "     units: %s" % i.units
        print 10*"-"

make_coord_system()