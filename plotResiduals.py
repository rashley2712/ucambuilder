#!/usr/bin/env python
import numpy
import matplotlib.pyplot
import matplotlib.image as mpimg
import argparse
import astropy.io.fits

def findClosestTime(data, target):
	distance = 1000.
	reading = (0, 0)
	for d in data:
		time = d[0]
		gap = time - target
		if abs(gap) < distance:
			distance = abs(gap)
			reading = (time, d[1])
	return reading
	
		

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Uses matplotlib to show the residuals of 2 light curves.')
	parser.add_argument('datafile1', type=str, help='Input data file for the light curves (usually a CSV file.')
	parser.add_argument('datafile2', type=str, help='Input data file for the light curves (usually a CSV file.')
	arg = parser.parse_args()
	print arg


	inputFile = open(arg.datafile1, 'r')
	
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
			except ValueError:
				print "No red value for MJD:", values[0]

			try:
				green = float(values[2])
				greens.append([time, green])
			except ValueError:
				print "No green value for MJD:", values[0]

			try:
				blue = float(values[3])
				blues.append([time, blue])
			except ValueError:
				print "No blue value for MJD:", values[0]

	# Discard the first and last reading
	del reds[0]
	del reds[-1]
	del greens[0]
	del greens[-1]	
	del blues[0]
	del blues[-1]	
	inputFile.close()
	
	rednp = numpy.array(reds)
	redvalues = rednp[:, 1]
	greennp = numpy.array(greens)
	greenvalues = greennp[:, 1]
	bluenp = numpy.array(blues)
	bluevalues = bluenp[:, 1]
	
	print redvalues
	print "Red mean:" , numpy.mean(redvalues)
	print "Red stddev:" , numpy.std(redvalues)
	print "Green mean:" , numpy.mean(greenvalues)
	print "Green stddev:" , numpy.std(greenvalues)
	print "Blue mean:" , numpy.mean(bluevalues)
	print "Blue stddev:" , numpy.std(bluevalues)

	inputFile = astropy.io.fits.open(arg.datafile2)

	alt_reds = [];
	headers = inputFile["CCD 1"].header
	data = inputFile["CCD 1"].data
	columns = inputFile["CCD 1"].columns
	objects = []
	
	for item in data:
		time = item[columns.names.index("MJD")]
		red = item[columns.names.index("Counts_2")]
		if (red!=0):
			alt_reds.append([time, red]);

	alt_greens = [];
	headers = inputFile["CCD 2"].header
	data = inputFile["CCD 2"].data
	columns = inputFile["CCD 2"].columns
	objects = []
	
	for item in data:
		time = item[columns.names.index("MJD")]
		green = item[columns.names.index("Counts_2")]
		if (green!=0):
			alt_greens.append([time, green]);

	alt_blues = [];
	headers = inputFile["CCD 3"].header
	data = inputFile["CCD 3"].data
	columns = inputFile["CCD 3"].columns
	objects = []
	
	for item in data:
		time = item[columns.names.index("MJD")]
		blue = item[columns.names.index("Counts_2")]
		if (blue!=0):
			alt_blues.append([time, blue]);
	
	# Discard the first and last reading
	del alt_reds[0]
	del alt_reds[-1]
	del alt_greens[0]
	del alt_greens[-1]	
	del alt_blues[0]
	del alt_blues[-1]	
	rednp = numpy.array(alt_reds)
	redvalues = rednp[:, 1]
	greennp = numpy.array(alt_greens)
	greenvalues = greennp[:, 1]
	bluenp = numpy.array(alt_blues)
	bluevalues = bluenp[:, 1]
	
	print redvalues
	print "Alt Red mean:" , numpy.mean(redvalues)
	print "Alt Red stddev:" , numpy.std(redvalues)
	print "Alt Green mean:" , numpy.mean(greenvalues)
	print "Alt Green stddev:" , numpy.std(greenvalues)
	print "Alt Blue mean:" , numpy.mean(bluevalues)
	print "Alt Blue stddev:" , numpy.std(bluevalues)
	
	
	MJDoffset = int(reds[0][0])
	
	print "MJD offset:",MJDoffset
	x_values = []
	y_values = []
	times = []
	for r in reds:
		x_values.append(r[0] - MJDoffset)
		times.append(r[0])
		y_values.append(r[1])

	for index, mjd in enumerate(times):
		(closestTime, counts) = findClosestTime(alt_reds, mjd)
		#print mjd, closestTime
		y_values[index] = y_values[index] - counts
		
	print "Res. Red mean:" , numpy.mean(y_values)
	print "Res. Red stddev:" , numpy.std(y_values)
		

	matplotlib.pyplot.plot(x_values, y_values, 'r.', label = 'i')
	matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset))
	matplotlib.pyplot.ylabel('Counts A - Counts T')

	x_values = []
	y_values = []
	times = []
	for g in greens:
		x_values.append(g[0] - MJDoffset)
		y_values.append(g[1])
		times.append(g[0])
		
	for index, mjd in enumerate(times):
		(closestTime, counts) = findClosestTime(alt_greens, mjd)
		#print mjd, closestTime
		y_values[index] = y_values[index] - counts
	print "Res. Green mean:" , numpy.mean(y_values)
	print "Res. Green stddev:" , numpy.std(y_values)
	
	matplotlib.pyplot.plot(x_values, y_values, 'g.', label = 'g')

	x_values = []
	y_values = []
	times = []
	for b in blues:
		x_values.append(b[0] - MJDoffset)
		y_values.append(b[1])
		times.append(b[0])
	
	for index, mjd in enumerate(times):
		(closestTime, counts) = findClosestTime(alt_blues, mjd)
		#print mjd, closestTime
		y_values[index] = y_values[index] - counts
	print "Res. Blue mean:" , numpy.mean(y_values)
	print "Res. Blue stddev:" , numpy.std(y_values)
	
		
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
	
	
	
