#!/usr/bin/env python


import argparse
import numpy as np
import methods
from numpy import linalg

parser = argparse.ArgumentParser(
	description = "Interface for spectrometer data analysis library")

parser.add_argument("filename",
	metavar	= "filename",
	type	= str,
	help	= "The filename of the original spectrometer data")

parser.add_argument("wavelength",
	metavar	= "nm",
	type	= int,
	help	= "The wavelength to display images at")

parser.add_argument("dimensions",
	metavar	= "pixel",
	type	= int,
	nargs	= 2,
	help	= "The dimensions of the image")

parser.add_argument("-o", "--original",
	action	= "store_true",
	help	= "Choice to display the original image")

parser.add_argument("-a", "--absorbance",
	action	= "store_true",
	help	= "Choice to display the aborbance image")

parser.add_argument("-s", "--svd",
	metavar	= "number of components",
	type	= int,
	default = 0,
	help	= "Required for -r, -u, -w, -v, and -d.")

parser.add_argument("-r", "--reverse",
	action	= "store_true",
	help	= """	Choice to display reverse svd image.
			Requires -s flag to be used.
			Required for -d.""")

parser.add_argument("-u",
	action	= "store_true",
	help	= """	Choice to display a U image for every component from svd.
			Requires -s flag to be used.""")

parser.add_argument("-w",
	action	= "store_true",
	help	= "Choice to display W matrix from svd. Requires -s flag to be used.")

parser.add_argument("-v",
	action	= "store_true",
	help	= """	Choice to display a Vt image for every component from svd.
			Requires -s flag to be used.""")

parser.add_argument("-d", "--difference",
	action	= "store_true",
	help	= """	Choice to display difference image between reverse svd image and
			original image. Requires -s and -r flags to be used.""")

parser.add_argument("-x", "--versus",
	metavar	= "component number",
	type	= int,
	nargs	= 2,
	default	= None,
	help	= "Takes in two component numbers to compare.")

args = parser.parse_args()

filename	= args.filename
wl		= args.wavelength
dimensions	= args.dimensions
original	= args.original
absorbance	= args.absorbance
components	= args.svd
reverse		= args.reverse
u		= args.u
w		= args.w
v		= args.v
difference	= args.difference
versus		= args.versus

x = dimensions[0]
y = dimensions[1]

data, nm = methods.get_data(filename)

data = data[:(x*y)]

if original:
  name = filename + "." + str(wl) + ".ori.eps"
  methods.show_image(methods.to2D(data, (x,y), wl), show=False, name=name, dir=filename)

if absorbance:
  name = filename + "." + str(wl) + ".abs.eps"
  absorbance = methods.absorbance(data, nm)
  methods.show_image(methods.to2D(absorbance, (x,y), wl), amp=True, show=False,
			name=name, dir=filename)

if components != 0:
  U, W, Vt = linalg.svd(data, full_matrices=0)

if reverse:
  name = filename +"." + str(wl) + ".rec.eps"
  reverse = methods.reverse_svd(U, W, Vt, components)
  methods.show_image(methods.to2D(reverse, (x,y), wl), show=False, name=name,
			dir=filename)
  
if u:
  first = True

  for i in range(components):
    name = "lsv" + str(i) + ".eps"
    if first:
      methods.show_image(methods.to2D(U, (x,y), i), show=False, name=name, dir=filename)
      first = False
    else:
      methods.show_image(methods.to2D(U, (x,y), i), amp=True, show=False,
			name=name, dir=filename)

if w:
  name = filename + ".sv.eps"
  methods.scatter_plot((range(1, len(W) + 1), W), xlog=True, ylog=True,
			highlight=components, left=0.9, show=False,
			name=name, dir=filename)

if v:
  first = True

  for i in range(components):
    name = "lsv" + str(i) + ".eps"
    if first:
      methods.show_image(methods.to2D(Vt, (y,x), i), show=False, name=name, dir=filename)
      first = False
    else:
      methods.show_image(methods.to2D(Vt, (y,x), i), amp=True, show=False,
			name=name, dir=filename)


if difference:
  name = filename + "." + str(wl) + ".dif.eps"
  diff = np.subtract(reverse, data)
  methods.show_image(methods.to2D(diff, (x,y), wl), amp=True, show=False,
			name=name, dir=filename)

#first arg: matrix, is either U or VT
#second arg: t, is the shared axis, should just be nm
#third arg: indexes, two element tuple containing the two component index you want to compare
#keyword arg: column, is True by default. Change to False if you want to compare rows of matricies for Vt

if versus is not None:
  name = filename + "." + str(versus[0]) + "." + str(versus[1]) + ".vs.vt.eps"
  methods.versus_graph(Vt, nm, (versus[0], versus[1]), column=False,
			name=name, dir=filename, show=False)

  name = filename + "." + str(versus[0]) + "." + str(versus[1]) + ".vs.u.eps"
  methods.versus_graph(U, nm, (versus[0], versus[1]), column=True,
			name=name, dir=filename, show=False)
