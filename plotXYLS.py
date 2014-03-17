#!/usr/bin/env python

import astropy.io.fits
import argparse
import matplotlib.pyplot, numpy
import Image, ImageDraw
import ultracamutils
import classes

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Takes a PNG file and an XYLS catalog as inputs and plots them on the screen or saves as a .png file')
	parser.add_argument('-i', '--input', type=str, help='Input file (XYLS)')
	parser.add_argument('-m', '--image', type=str, help='Input file (PNG)')
	parser.add_argument('-o', '--output', default='testplot.png', type=str, help='Output file (PNG)')
	parser.add_argument('-p', '--preview', action='store_true', help='Show a preview on screeen using Matplotlib')
	parser.add_argument('-d', '--debuglevel', default = 1, type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	arg = parser.parse_args()

	if (arg.input!=None): 
		inputFilename = arg.input
	else:
		print "Please specify an input file (XYLS catalog)"
		exit()

	if (arg.image!=None): 
		imageFilename = arg.image
	else:
		print "Please specify an input file (image)"
		exit()

	outputFilename = arg.output
	
	debug = classes.debugObject(arg.debuglevel)
	
	""" Load the PNG image first
	"""
	#imageHDUList = astropy.io.fits.open(fitsFilename)
	#imageData = imageHDUList[0].data
	#imageHDUList.close()

	imageData=matplotlib.image.imread(imageFilename)

	#imageData = ultracamutils.percentiles(imageData, 50, 90)
	imageData = numpy.flipud(imageData)
	imgplot = matplotlib.pyplot.imshow(imageData, cmap="gray")

	""" Load the catalog file
	"""
	tableHDUList = astropy.io.fits.open(inputFilename)
	tableData = tableHDUList[1].data
	tableColumns = tableHDUList[1].columns
	catalog = []
	for item in tableData:
		object = {'id': 0, 'x': 0, 'y': 0, 'flux': 0}
		object['id'] = item[tableColumns.names.index("ID")]
		object['x'] = item[tableColumns.names.index("X")]
		object['y'] = item[tableColumns.names.index("Y")]
		object['flux'] = item[tableColumns.names.index("FLUX")]
		catalog.append(object)
		
	print catalog
	for i in catalog:
		x = i['x'] 
		y = i['y'] 
		fwhm = 4
		matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), fwhm*3, color='green', fill=False, linewidth=2.0))

	matplotlib.pyplot.show()
