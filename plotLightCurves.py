#!/usr/bin/env python
import numpy
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Uses matplotlib to make nice light curves.')
	parser.add_argument('datafile', type=str, help='Input data file for the light curves (usually a CSV file.')
	arg = parser.parse_args()
	print arg

	inputFile = open(arg.datafile, 'r')
	
	times = [];
	reds = [];
	
	for n, line in enumerate(inputFile):
		values = line.split(',')
		if n==0:
			headings = values
		else:
			
			time = float(values[0])			
			try:
				red = float(values[1])
				reds.append([time, red])
			except ValueError:
				print "No red value for MJD:", values[0]
			
	inputFile.close()

	MJDoffset = int(reds[0][0])
	print "MJD offset:",MJDoffset
	x_values = []
	y_values = []
	for r in reds:
		x_values.append(r[0] -MJDoffset)
		y_values.append(r[1])

	
	matplotlib.pyplot.plot(x_values, y_values, 'r.')
	#c='r', marker=',', linestyle='none')
	matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset))
	matplotlib.pyplot.ylabel('Counts')
	
	
	# Now check everything with the defaults:
	fig = matplotlib.pyplot.gcf()
	DPI = fig.get_dpi()
	print "DPI:", DPI
	DefaultSize = fig.get_size_inches()
	print "Default size in Inches", DefaultSize
	print "Which should result in a %i x %i Image"%(DPI*DefaultSize[0], DPI*DefaultSize[1])

	matplotlib.pyplot.show()
	
	fig.savefig('test2.eps',dpi=100, format='eps')
	
	
	
