#!/usr/bin/env python

import astropy.io.fits
import argparse
import matplotlib.pyplot, numpy, math
import matplotlib.image
import Image, ImageDraw
import ultracamutils
import classes, wcsclasses
import os, subprocess, sys, json

def getArrayFromObjects(objects, propertyName):
	values = []
	for o in objects:
		value = o[propertyName]
		values.append(value)
	return numpy.array(values)

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Compares the wcs solutions produced by Astrometry.net across r, g and b.')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-d', '--debuglevel', default = 1, type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	
	arg = parser.parse_args()

	config = ultracamutils.readConfigFile(arg.configfile)
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);

	""" First check if a directory is made for the output files in the working director folder. and create one. 
	"""
	runDate, runNumber = ultracamutils.separateRunNameAndDate(arg.runname)
	
	channels = ['r', 'g', 'b']
	channelDescriptions = {'r':'Red', 'g':'Green', 'b':'Blue'}
	
		
	""" Read info about the run to get the starting RA and DEC coordinates
	"""
	debug.write("Getting run info from the file:" + config.RUNINFO, level = 2)
	runInfo = ultracamutils.getRunInfo(config.RUNINFO, arg.runname)
	
	debug.write("Run Info:\n----------------------", level = 2)
	debug.write(runInfo, level = 2)
	debug.write("----------------------", level = 2)

	solved = [False, False, False]
	for n,c in enumerate(channels):
		# Check that there is a solution for each channel
		solvedFileMarker = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + '_' + c + ".solved"
		if os.path.exists(solvedFileMarker):
			solved[n] = True
	print solved
	if solved[0] & solved[1] & solved[2]:
		print "All solved!"

	
	""" Load the fits files of the wcs solutions
	"""
	wcsSolutions = {'r': None, 'g': None, 'b': None }
	for n,c in enumerate(channels):
		if solved[n]:
			wcsFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + '_' + c + ".wcs"
			wcsJSONFilename = ultracamutils.addPaths(config.SITE_PATH, arg.runname) + '_' + c + "_wcs.json"
			wcsFile = astropy.io.fits.open(wcsFilename)
			header = wcsFile[0].header
			wcs = wcsclasses.wcsSolution()
	
			equinox = float(header['EQUINOX'])
			referenceCoord = float(header['CRVAL1']), float(header['CRVAL2'])
			referencePixel = float(header['CRPIX1']), float(header['CRPIX2'])
			
			CD_array = [ [header['CD1_1'], header['CD1_2']], [ header['CD2_1'], header['CD2_2'] ] ]
				
			wcs.setSolution(equinox, referenceCoord, referencePixel, CD_array)
			
			aOrder = int(header['A_ORDER'])
			for i in range(aOrder+1): 
				for j in range(aOrder+1):
					if (i+j)>aOrder+1: continue;
					paramString = "A_" + str(i) + "_" + str(j)
					sip = paramString
					try:
						sipvalue = float(header[sip])
						wcs.setSIPAvalue(sip, i, j, sipvalue)
					except KeyError:
						sipvalue = 0.0

			bOrder = int(header['B_ORDER'])
			for i in range(bOrder+1): 
				for j in range(bOrder+1):
					if (i+j)>bOrder+1: continue;
					paramString = "B_" + str(i) + "_" + str(j)
					sip = paramString
					try:
						sipvalue = float(header[sip])
						wcs.setSIPBvalue(sip, i, j, sipvalue)
					except KeyError:
						sipvalue = 0.0

			aOrder = int(header['AP_ORDER'])
			for i in range(aOrder+1): 
				for j in range(aOrder+1):
					if (i+j)>aOrder+1: continue;
					paramString = "AP_" + str(i) + "_" + str(j)
					sip = paramString
					try:
						sipvalue = float(header[sip])
						wcs.setSIPAPvalue(sip, i, j, sipvalue)
					except KeyError:
						sipvalue = 0.0

			bOrder = int(header['BP_ORDER'])
			for i in range(bOrder+1): 
				for j in range(bOrder+1):
					if (i+j)>bOrder+1: continue;
					paramString = "BP_" + str(i) + "_" + str(j)
					sip = paramString
					try:
						sipvalue = float(header[sip])
						wcs.setSIPBPvalue(sip, i, j, sipvalue)
					except KeyError:
						sipvalue = 0.0
					

			print channelDescriptions[c]
			print wcs
			wcsSolutions[c] = wcs
			
	
	# Test case
	print "Test: Going from pixel to world"
	pixel = (466.39, 482.31)
	world = wcsSolutions['r'].getWorldSIP(pixel)
	print "Test pixel:", pixel
	print "Test world:", world
	print
	
	print "Test: Going from world to pixel"
	world = (296.04, 40.29)
	pixel = wcsSolutions['r'].getPixel(world)
	print "Test world:", world
	print "Test pixel:", pixel
	
	colour = 'b'
	catalogFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_" + colour + ".xyls"
	print "Now load the red input catalog:", catalogFilename
	catalogFile = astropy.io.fits.open(catalogFilename)
	headers = catalogFile[1].header
	columns = catalogFile[1].columns
	data = catalogFile[1].data
	
	print catalogFile.info()
	objects = []
	for d in data:
		object = {}
		ID = d[columns.names.index('ID')]
		x = d[columns.names.index('X')]
		y = d[columns.names.index('Y')]
		flux = d[columns.names.index('FLUX')]
		print ID, x, y, flux
		object['ID'] = ID
		object['x'] = x
		object['y'] = y
		object['flux'] = flux
		objects.append(object)

	figure = matplotlib.pyplot.figure(figsize=(12, 12))
	
	ymax = 1032
	arrowScale = 10000000
	for o in objects:
		# First transform to world coordinates without using SIP
		(x, y) = o['x'], o['y']
		(world_x, world_y) = wcsSolutions['r'].getWorldCoord((x,y))
		print x, y, ":", world_x, world_y, 
		# Now revert back to pixel position with SIP polynomial
		(new_x, new_y) = wcsSolutions['r'].getPixel((world_x, world_y))
		x_arrow = new_x - x
		y_arrow = new_y - y
		print x_arrow, y_arrow
		matplotlib.pyplot.plot([x, x + x_arrow * arrowScale], [ ymax - y, ymax - (y + y_arrow * arrowScale)], lw=1, color='black')
	
	
	matplotlib.pyplot.gca().invert_yaxis()
	
	# Also load the image bitmap
	pngFile = ultracamutils.addPaths(config.SITE_PATH, arg.runname) + '_' + colour + ".png"
	print "Loading bitmap:", pngFile
	img = matplotlib.image.imread(pngFile)
	
	imgplot = matplotlib.pyplot.imshow(img)
	
	matplotlib.pyplot.show()
	
	figure.savefig('test2.eps',dpi=100, format='eps')
	
	"""
	Now compare the offsets from colour to colour
	"""
	
	"""
	gridSize = 1024
	gridSpacing = 20
	arrowScale = 2.0
	greenOffsetx = numpy.zeros((gridSize, gridSize))
	greenOffsety = numpy.zeros((gridSize, gridSize))
	blueOffsetx = numpy.zeros((gridSize, gridSize))
	blueOffsety = numpy.zeros((gridSize, gridSize))
	
	for i in range(0, gridSize, gridSpacing):
		print i
		for j in range(0, gridSize, gridSpacing): 
			greenPixel = (i, j)
			world = wcsSolutions['g'].getWorldSIP(greenPixel)
			redPixel = wcsSolutions['r'].getPixel(world)
			offset = (greenPixel[0] - redPixel[0], greenPixel[1] - redPixel[1])
			greenOffsetx[i][j] = offset[0]
			greenOffsety[i][j] = offset[1]
			matplotlib.pyplot.plot([i, i+offset[0] * arrowScale], [j, j+offset[1] * arrowScale], lw=1, color='green')
			
			bluePixel = (i, j)
			world = wcsSolutions['b'].getWorldSIP(greenPixel)
			redPixel = wcsSolutions['r'].getPixel(world)
			offset = (bluePixel[0] - redPixel[0], bluePixel[1] - redPixel[1])
			blueOffsetx[i][j] = offset[0]
			blueOffsety[i][j] = offset[1]
			matplotlib.pyplot.plot([i, i + offset[0] * arrowScale ], [j, j + offset[1] * arrowScale ], lw=1, color='blue')
			
		
	print greenOffsetx
	print greenOffsety
			
	print blueOffsetx
	print blueOffsety

	matplotlib.pyplot.show()
	"""

	"""X,Y = numpy.meshgrid( numpy.arange(0,gridSize,1), numpy.arange(0, gridSize, 1) )
	U = numpy.cos(X)
	V = numpy.sin(Y)
	
	#QP = matplotlib.pyplot.quiver(X, Y, greenOffsetx, greenOffsety, scale = 1.5)
	#QP = matplotlib.pyplot.quiver(X, Y, U, V, scale=2)
	#matplotlib.pyplot.quiverkey(QP, """