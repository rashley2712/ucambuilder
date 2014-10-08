#!/usr/bin/env python
import numpy
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse
import astropy.io.fits

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Uses matplotlib to make nice light curves.')
	parser.add_argument('datafile', type=str, help='Input data file for the light curves (usually a CSV file.')
	parser.add_argument('-f', action='store_true', help='Input file is a FITS file, not a CSV.')
	arg = parser.parse_args()
	print arg


	if (not arg.f):
		inputFile = open(arg.datafile, 'r')
		
		times = [];
		reds = [];
		greens = [];
		blues = [];
		
		for n, line in enumerate(inputFile):
			values = line.split(',')
			if n==0:
				headings = values
			else:
				
				time = float(values[0])			
				try:
					red = float(values[1])
					reds.append([time, red])
				except:
					print "No red value for MJD:", values[0]
	
				try:
					green = float(values[2])
					greens.append([time, green])
				except:
					print "No green value for MJD:", values[0]
	
				try:
					blue = float(values[3])
					blues.append([time, blue])
				except:
					print "No blue value for MJD:", values[0]
	
				
		inputFile.close()
	else:
		inputFile = astropy.io.fits.open(arg.datafile)
	
		reds = [];
		headers = inputFile["CCD 1"].header
		data = inputFile["CCD 1"].data
		columns = inputFile["CCD 1"].columns
		objects = []
		
		for item in data:
			time = item[columns.names.index("MJD")]
			red = item[columns.names.index("Counts_2")]
			if (red!=0):
				reds.append([time, red]);
	
		greens = [];
		headers = inputFile["CCD 2"].header
		data = inputFile["CCD 2"].data
		columns = inputFile["CCD 2"].columns
		objects = []
		
		for item in data:
			time = item[columns.names.index("MJD")]
			green = item[columns.names.index("Counts_2")]
			if (green!=0):
				greens.append([time, green]);

		blues = [];
		headers = inputFile["CCD 3"].header
		data = inputFile["CCD 3"].data
		columns = inputFile["CCD 3"].columns
		objects = []
		
		for item in data:
			time = item[columns.names.index("MJD")]
			blue = item[columns.names.index("Counts_2")]
			if (blue!=0):
				blues.append([time, blue]);
	
			
	
	MJDoffset = int(reds[0][0])
	
	print "MJD offset:",MJDoffset
	x_values = []
	y_values = []
	for r in reds:
		x_values.append(r[0] -MJDoffset)
		y_values.append(r[1])

	matplotlib.pyplot.plot(x_values, y_values, 'r.', label = 'i')
	matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset))
	matplotlib.pyplot.ylabel('Counts')

	x_values = []
	y_values = []
	for g in greens:
		x_values.append(g[0] -MJDoffset)
		y_values.append(g[1])
	matplotlib.pyplot.plot(x_values, y_values, 'g.', label = 'g')

	x_values = []
	y_values = []
	for b in blues:
		x_values.append(b[0] -MJDoffset)
		y_values.append(b[1])
	matplotlib.pyplot.plot(x_values, y_values, 'b.', label = 'u')

	matplotlib.pyplot.legend()
	
	# Now check everything with the defaults:
	fig = matplotlib.pyplot.gcf()
	DPI = fig.get_dpi()
	print "DPI:", DPI
	DefaultSize = fig.get_size_inches()
	print "Default size in Inches", DefaultSize
	print "Which should result in a %i x %i Image"%(DPI*DefaultSize[0], DPI*DefaultSize[1])

	matplotlib.pyplot.show()
	
	fig.savefig('test2.eps',dpi=100, format='eps')
	
	
	
