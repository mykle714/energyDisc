#!/usr/bin/env python

usage = """

"""

import os, sys, time, tempfile, argparse
from math    import *
from Gnuplot import Gnuplot
from cpl     import mathematics, statistics
from cpl.mac import vm, amoeba
from cpl.ccl import symm, exp, det, utility

def bindex(idx, bin):
   """Given a 0-offset index idx of a binXbin binned pixel, returns the
      center coordinate of the binned pixel.  The returned coordinate is with
      respect to a system in which each unbinned pixel centers on its 1-offset
      sequence number.  For example, the center of the first pixel is at 1,
      and the center of the n-th pixel is at n.
   """
   return idx*bin + 0.5*bin + 0.5

def view(FORMAT, image, Bin        = 1,
                        Zoom       = 100, # if 0, window 1 only
                        Range      = (0, 0),
                        color      = False,
                        logarithm  = True,
                        Center     = (0, 0),
                        distance   = 0,
                        pixel      = 0,
                        wavelength = (0, 0),
                        ell        = '',
                        spt        = '',
                        symbolsize = 1,
                        eps        = '',
                        epssize    = 7,
                        fontsize   = 12):
   intmin, intmax = Range
   lmin, lmax     = wavelength
   auto     = intmin*intmax == 0
   sample   = []
   Image    = image.getTitle()
   sizex    = image.getSizeX()
   sizey    = image.getSizeY()
   sizez    = image.saturate()
   binned   = tempfile.mktemp(".dat")
   data     = open(binned, 'w')
   for i in range(sizex / Bin):
      x = bindex(i, Bin)
      for j in range(sizey / Bin):
         y = bindex(j, Bin)
         z = image.intensity(x, y, Bin)
         if auto and z > 1 and sizex/30 < i < sizex/3:
            sample.append(2 * int(0.5*z))
         data.write("%g %g %g\n" % (x, y, z))
      data.write("\n")
   data.close()
   intmin = intmin if intmin else mathematics.tickmax(min(sample))
   intmax = intmax if intmax else mathematics.tickmax(
                                     2*statistics.mode(sample) + intmin)
   if intmin >= intmax: intmax = intmin + 1
   zoomin = tempfile.mktemp(".dat")

   gp = Gnuplot(persist = False)
   gp("unset key")
   gp("set   cblabel 'Detector count'")
   gp("set   view map")
   gp("set   size ratio -1")
   gp("set   palette %s" % ("color" if color else "gray negative"))
   gp("set   cbrange [%g:%g]" % (intmin, intmax))
   if logarithm: gp("set logscale cb")
   gp("set   grid front lc rgb 'royalblue'")
   xlabel = "Click inside to zoom in or outside to exit." if Zoom else\
            "Detector pixel"
   window1(gp, Image, binned, [], 0, sizex, sizey, Bin, intmin, intmax,
           Center[0], Center[1], distance, pixel, lmin, lmax, ell, spt,
           symbolsize, xlabel, eps, epssize, fontsize)

   while Zoom:
      if os.path.isfile(zoomin): os.remove(zoomin)
      gp("pause mouse button1")
      gp("u  = MOUSE_X")
      gp("v  = MOUSE_Y")
      gp("set print '%s'" % zoomin)
      gp("print u,v")
      gp("set print")
      while not os.path.isfile(zoomin) or not os.access(zoomin, os.R_OK):
         time.sleep(0.1)
      data = open(zoomin)
      uv   = data.readline().split()
      u, v = int(float(uv[0]) + 0.5), int(float(uv[1]) + 0.5)
      data.close()
      if image.outside(u, v): sys.exit()

      window2(gp, Image, image,  (u, v), Zoom, sizex, sizey, Bin,
              Center[0], Center[1], distance, pixel, lmin, lmax, ell, spt,
              symbolsize)
      window1(gp, Image, binned, (u, v), Zoom, sizex, sizey, Bin,
              intmin, intmax, Center[0], Center[1], distance, pixel,
              lmin, lmax, ell, spt, symbolsize, xlabel)

   return gp, binned, (intmin, intmax)

def window1(gp, Image, binned, zoomat, zmw, sizex, sizey, bin, intmin, intmax,
            centerx, centery, distance = 0, pixel = 0, lmin = 0, lmax = 0,
            ell = '', spt = '', symbolsize = 1, xlabel = '', eps = '',
            epssize = 7, fontsize = 12):
   if not gp: return
   gp("set terminal x11 1 title '%s' size %d, %d enhanced font 'Helvetica,%d'"
      % (Image, 1.15 * sizex/bin, sizey/bin, 12))
   gp("set   colorbox")
   if not (log10(intmin)).is_integer():
      gp("set   cbtics add ('%g' %g)" % (intmin, intmin))
   if not (log10(intmax)).is_integer():
      gp("set   cbtics add ('%g' %g)" % (intmax, intmax))
   gp("x0 = %g" % centerx)
   gp("y0 = %g" % centery)
   gp("ps = %g" % pixel)
   gp("set xrange [1:%d]" % sizex)
   gp("set yrange [%d:1]" % sizey)
   gp("set trange [-pi:pi]")
   if xlabel: gp("set xlabel '%s'" % xlabel)
   if zmw:
      u, v = zoomat
      gp("set object 1 rect at first %g,%g size first %g,%g front fs empty border rgb 'royalblue'"
         % (u, v, zmw, zmw))
   if abs(complex(centerx, centery)) and (lmin or lmax):
      obj = 2
      for d in (5, 3, 2, 1.5, 1.2):
         rmin = distance * tan(2 * asin(0.5 * lmin / d)) / pixel
         rmax = distance * tan(2 * asin(0.5 * lmax / d)) / pixel
         gp("set object %d circle at first x0,y0 size first %g front fs empty border rgb 'royalblue'"
            % (obj, rmin))
         obj += 1
         if lmin != lmax:
            gp("set object %d circle at first x0,y0 size first %g front fs empty border rgb 'light-red'"
               % (obj, rmax))
            obj += 1
   gp("set parametric")
   cmd  = "plot "
   cmd += "'%s' u 1:2:3 w image, " % binned
   if abs(complex(centerx, centery)):
      cmd += "x0,y0 w p pt 1 ps 2 lc rgb 'royalblue', "
      if ell and distance and pixel:
         data = open(ell)
         for i in data:
            line = i.split()
            x    = (float(line[3]) - centerx) * pixel 
            y    = (float(line[4]) - centery) * pixel
            a    = sqrt(x*x + y*y)
            e    = (-distance + sqrt(distance*distance + 4*a*a)) / (2*a)
            beta = sqrt(1-e*e)
            cmd += "(%g*(1+cos(t))-%g*sin(t)*%g)/ps+x0, " % (x, y, beta)
            cmd +=\
   "(%g*(1+cos(t))+%g*sin(t)*%g)/ps+y0 w l lt 1 lw 1 lc rgb 'light-red', " %\
   (y, x, beta)
         data.close()
   if spt:
      cmd += "'%s' u 4:5 w p pt 2 ps %g lc rgb 'light-red', " % (spt,
                                                                 symbolsize)
   gp(cmd[:-2])

   if eps:
      gp("set output '%s'" % eps)
      gp("set encoding iso_8859_1")
      gp(
   "set terminal postscript eps enhanced color solid 'Helvetica' %d size %f,%f"
   % (fontsize, 1.15*epssize, epssize))
      gp("unset xlabel")
      gp("replot")
      gp("set output")
      gp(
   "set terminal x11 1 title '%s' size %d, %d enhanced font 'Helvetica,%d'"
   % (Image, 1.15 * sizex/bin, sizey/bin, 12))
   gp("unset parametric")

def window2(gp, Image, image, zoomat, zmw, sizex, sizey, bin,
            centerx, centery, distance = 0, pixel = 0, lmin = 0, lmax = 0,
            ell = '', spt = '', symbolsize = 1):
   if not gp: return
   u, v = zoomat
   size = min(sizex, sizey) / bin
   gp("set terminal x11 2 title '%s' size %d, %d enhanced font 'Helvetica,%d'"
      % (Image, size, size, 12))
   gp("unset colorbox")
   gp("unset xlabel")
   gp("unset object")
   gp("x0 = %g" % centerx)
   gp("y0 = %g" % centery)
   gp("ps = %g" % pixel)
   gp("set xrange [%d:%d]" % (u-zmw/2, u+zmw/2))
   gp("set yrange [%d:%d]" % (v+zmw/2, v-zmw/2))
   gp("set parametric")
   cmd  = "plot "
   cmd += "'-' u 1:2:3 w image, "
   if abs(complex(centerx, centery)):
      cmd += "x0,y0 w p pt 1 ps 2 lc rgb 'royalblue', "
      if ell and distance and pixel:
         data = open(ell)
         for i in data:
            line = i.split()
            x    = (float(line[3]) - centerx) * pixel
            y    = (float(line[4]) - centery) * pixel
            a    = sqrt(x*x + y*y)
            e    = (-distance + sqrt(distance*distance + 4*a*a)) / (2*a)
            beta = sqrt(1-e*e)
            cmd += "(%g*(1+cos(t))-%g*sin(t)*%g)/ps+x0, " % (x, y, beta)
            cmd +=\
   "(%g*(1+cos(t))+%g*sin(t)*%g)/ps+y0 w l lt 1 lw 1 lc rgb 'light-red', " %\
   (y, x, beta)
         data.close()
   if spt:
      cmd += "'%s' u 4:5 w p pt 2 ps %g lc rgb 'light-red', " % (spt,
                                                                 symbolsize)
   gp(cmd[:-2])
   for i in range(u - zmw/2, u + zmw/2 + 1):
      for j in range(v + zmw/2, v - zmw/2 - 1, -1):
         gp("%d %d %g" % (i, j, image.intensity(i, j)
                             if image.inside(i, j) else 0))
      gp('')
   gp('e')
   gp("unset parametric")

if __name__ == "__main__":
   Format = sys.argv[0].split('/')[-1][:-3]
   defbin = [2, 100]
   if Format == "RayonixMX340": defbin = [6, 300]
   if Format == "MARCCD165":    defbin = [4, 200]
   if Format == "IMAGINE":      defbin = [5, 100]

   parser = argparse.ArgumentParser(
      description = "Display a diffraction image in %s format." % Format)
   parser.add_argument("filename",
      metavar = "filename",
      type    = str,
      help    = "The filename of the diffraction image to be displayed")
   parser.add_argument("-b", "--binning",
      metavar = "pixel",
      type    = int,
      nargs   = 2,
      default = defbin,
      help    = """The entire image is NxN binned and displayed in the first
                   window, where N is given by the first number.  A part of MxM
                   pixels in the image is displayed in the second window, where
                   M is given by the second number.""")
   parser.add_argument("-d", "--detector",
      metavar = "data",
      type    = float,
      nargs   = 4,
      default = [0, 0, 0, 0],
      help    = """Detector parameters: direct-beam-center x, y in pixel,
                   crystal-to-detector distance in mm, detector pixel size in
                   mm""")
   parser.add_argument("-w", "--wavelength",
      metavar = "lambda",
      type    = float,
      nargs   = '+',
      default = [0, 0],
      help    = """Wavelength or wavelength range for resolution markings.
                   Detector parameters are also required to generate resolution
                   markings at 5, 3, 2, 1.5, and 1.2 Angstrom.  Blue and red
                   circles correspond to the specified short and long
                   wavelengths, respectively.""")
   parser.add_argument("-c", "--color",
      action  = "store_true",
      help    = "Use false colors.")
   parser.add_argument("-g", "--logarithmic",
      action  = "store_true",
      help    = "Display intensity in logrithmic scale.")
   parser.add_argument("-r", "--range",
      metavar = "intensity",
      type    = int,
      nargs   = 2,
      default = [0, 0],
      help    = "The range of detector intensity to be displayed")
   parser.add_argument("-e", "--ellipses",
      metavar = "sptFile",
      type    = str,
      default = '',
      help    = "A .spt filename to overlay Laue ellipses")
   parser.add_argument("-s", "--spots",
      metavar = "sptFile",
      type    = str,
      default = '',
      help    = "A .spt filename to overlay spots")
   parser.add_argument("-y", "--symbolsize",
      metavar = "pointsize",
      type    = float,
      default = 1,
      help    = "The symbol size for marking the spots")
   parser.add_argument("-f", "--fontsize",
      metavar = "point",
      type    = int,
      default = 12,
      help    = "The font size for encapsulated PostScript output")
   parser.add_argument("-l", "--length",
      metavar = "inch",
      type    = float,
      default = 7,
      help    = "The length in inch for encapsulated PostScript output")
   parser.add_argument("-o", "--output",
      metavar = "eps",
      type    = str,
      default = '',
      help    = "The filename for encapsulated PostScript output")
   args = parser.parse_args()

   Image      = args.filename
   bin, zmw   = args.binning
   ctr        = args.detector[:2]
   distance   = args.detector[2]
   pixel      = args.detector[3]
   lmin       = min(args.wavelength[:2])
   lmax       = max(args.wavelength[:2])
   color      = args.color
   logarithm  = args.logarithmic
   ran        = args.range
   ell        = args.ellipses
   spt        = args.spots
   symbolsize = args.symbolsize
   fontsize   = args.fontsize
   epssize    = args.length
   eps        = args.output

   FORMAT     = eval("exp." + Format)
   image      = FORMAT(Image)
   sizex      = image.getSizeX()
   sizey      = image.getSizeY()
   while sizex%bin or sizey%bin: bin -= 1

   view(FORMAT, image, bin, zmw, ran, color, logarithm, ctr, distance, pixel,
        (lmin, lmax), ell, spt, symbolsize, eps, epssize, fontsize)
