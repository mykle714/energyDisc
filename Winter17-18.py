import sys
import numpy as np
import methods, time
from numpy import linalg

#filename = sys.argv[1]
filename = "Hamilton10_40um"
wavelength = 900 #nm
cutoff = 40000

start_time = int(time.time())

data, nm = methods.get_data(filename, cutoff)

mid_time = int(time.time())

print("read time: %d" % (mid_time - start_time))

wl = methods.get_wavelength(nm, wavelength)

name = filename + "." + str(wavelength) + ".ori.eps"

dark_aligned = methods.align_dark(data)

single = methods.isolate_wavelength(dark_aligned, wl)

sqared = methods.square(single)

del(data)

methods.show_image(sqared)

end_time = int(time.time())

print("total time: %d" % (end_time - start_time))