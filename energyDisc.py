import sys
import numpy as np
import methods
from numpy import linalg

#filename = sys.argv[1]
filename = "Washington"
wl = 1000
comp = 2

data, nm = methods.get_data(filename)       #your code

name = filename + "." + str(wl) + ".ori.eps"

methods.show_image(methods.to2D(data, (41,40), wl), show=False, name=name, dir=filename)         #dislays the image at a particular wavelength

#to show multiple figures at once, set keyword arg, show, to False except the last figure you would like to sho
#show is True by default
#methods.show_image(methods.to2D(data1, (41,40), wl), show=False)
#methods.show_image(methods.to2D(data2, (41,40), wl))

#add two arguments: name and subdirectory to save. All subdirectories and images will appear in the ./figures/ directory
#methods.show_image(methods.to2D(data, (41,40), wl), "washington_img", "./imgs")

name = filename + "." + str(wl) + ".abs.eps"

absorbance = methods.absorbance(data, nm)       #gets the (anti-)aborbance ratio from the data
methods.show_image(methods.to2D(absorbance, (41,40), wl), amp=True, show=False, name=name, dir=filename)

U, W, VT = linalg.svd(data, full_matrices=0)     #svd
first = True

for i in range(comp):
    name = "lsv" + str(i) + ".eps"
    if first:
        methods.show_image(methods.to2D(U, (41,40), i), show=False, name=name, dir=filename)
        first = False
    else:
        methods.show_image(methods.to2D(U, (41,40), i), amp=True, show=False, name=name, dir=filename)


#isolates the first 10 components and returns the data in the original format (the format that methods.get_data() returns)

name = filename +"." + str(wl) + ".rec.eps"
reverse = methods.reverse_svd(U, W, VT, comp)
methods.show_image(methods.to2D(reverse, (41,40), wl), show=False, name=name, dir=filename)

name = filename + ".sv.eps"
methods.scatter_plot((range(1, len(W) +1), W), xlog=True, ylog=True, highlight=comp, left=0.9, show=False, name=name, dir=filename)

name = filename + "." + str(wl) + ".dif.eps"
diff = np.subtract(reverse, data)
methods.show_image(methods.to2D(diff, (41,40), wl), amp=True, show=False, name=name, dir=filename)

#first arg: matrix, is either U or VT
#second arg: t, is the shared axis, should just be nm
#third arg: indexes, two element tuple containing the two component index you want to compare
#keyword arg: column, is True by default. Change to False if you want to compare rows of matricies for Vt
x = 0
y = 1
name = filename + "." + str(x) + "." + str(y) + ".vs.eps"
methods.versus_graph(VT, nm, (x,y), column=False, name=name, dir=filename, show=False)