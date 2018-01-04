#!/usr/bin/env python

import methods
import argparse

parser = argparse.ArgumentParser(
	description = "Diagnostics for spectrometer data.")

parser.add_argument("filename",
	metavar	= "filename",
	type	= str,
	help	= "filename")

args = parser.parse_args()

filename = args.filename

data, nm = methods.get_data(filename)

print "length of data: " + str(len(data))
