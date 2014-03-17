#!/usr/bin/env python
import ultracamutils
import sys, subprocess, re, json
import classes, numpy as np
import astropy.io.fits
import argparse
	
if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(description='Reads the files produced by earlier steps in the pipeline.')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	parser.add_argument('--xyls', action='store_true', help='Create an XY-list for each channel. Used by Astrometry.net')
	parser.add_argument('-d', '--debuglevel', type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	arg = parser.parse_args()
	
	config = ultracamutils.readConfigFile(arg.configfile)
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);
	
	runName = sys.argv[1]
	
	redFilename = ultracamutils.addPaths(config.SITE_PATH, runName) + "_r.json"	
	redObjects = ultracamutils.buildObjectsFromJSON(redFilename)
	
	greenFilename = ultracamutils.addPaths(config.SITE_PATH, runName) + "_g.json"
	greenObjects = ultracamutils.buildObjectsFromJSON(greenFilename)
	
	blueFilename = ultracamutils.addPaths(config.SITE_PATH, runName) + "_b.json"
	blueObjects = ultracamutils.buildObjectsFromJSON(blueFilename)
	
	totalRedObjects = len(redObjects)
	totalGreenObjects = len(greenObjects)
	totalBlueObjects = len(blueObjects)
		
	redObjects = ultracamutils.filterOutCosmicRays(redObjects)
	greenObjects = ultracamutils.filterOutCosmicRays(greenObjects)
	blueObjects = ultracamutils.filterOutCosmicRays(blueObjects)

	redObjects = ultracamutils.filterOutLowFrameCountObjects(redObjects, 50)
	greenObjects = ultracamutils.filterOutLowFrameCountObjects(greenObjects, 50)
	blueObjects = ultracamutils.filterOutLowFrameCountObjects(blueObjects, 50)
	
	debug.write("%d red objects after cosmic ray and low frames filtering, was %d"%(len(redObjects), totalRedObjects), level = 2)
	debug.write("%d green objects after cosmic ray and low frames filtering, was %d"%(len(greenObjects), totalGreenObjects), level = 2)
	debug.write("%d blue objects after cosmic ray and low frames filtering, was %d"%(len(blueObjects), totalBlueObjects), level = 2)

	print "Timing check red"
	redObjects = ultracamutils.filterOutBadTimingFrames(redObjects)
	print "Timing check green"
	greenObjects = ultracamutils.filterOutBadTimingFrames(greenObjects)
	print "Timing check blue"
	blueObjects = ultracamutils.filterOutBadTimingFrames(blueObjects)
	
	
		
	for i in redObjects:
		meanPosition = i.calculateMeanPosition()
		meanFlux = i.calculateMeanFlux()
	for i in greenObjects:
		meanPosition = i.calculateMeanPosition()
		meanFlux = i.calculateMeanFlux()
	for i in blueObjects:
		meanPosition = i.calculateMeanPosition()
		meanFlux = i.calculateMeanFlux()
	
	if (arg.xyls):
		channelNames = ['r','g','b']
		debug.write("Creating an XYLS file")
		
		sortedObjects = sorted(blueObjects, key= lambda p: p.meanFlux, reverse=True)
		x_values, IDs, y_values, fluxes = [], [], [], []
		
		totalObjects = len(sortedObjects)
		if totalObjects<10:
			outputLength = totalObjects
		elif totalObjects<20:
			outputLength = totalObjects * 0.8
		elif totalObjects<50:
			outputLength = totalObjects * 0.7
		elif totalObjects<100:
			outputLength = totalObjects * 0.5
		else:
			outputLength = totalObjects * 0.1
		outputLength = int(outputLength)
		
		for num, i in enumerate(sortedObjects):
			IDs.append(i.id)
			x_values.append(i.meanPosition[0])
			y_values.append(i.meanPosition[1])
			fluxes.append(i.meanFlux)
			if num==outputLength: break;

		print fluxes

		FITSFilename = "test_blue.xyls"
		debug.write("Writing FITS file: " + FITSFilename, level=1)
	
		col1 = astropy.io.fits.Column(name='ID', format='I', array=IDs)
		col2 = astropy.io.fits.Column(name='X', format='E', array=x_values)
		col3 = astropy.io.fits.Column(name='Y', format='E', array=y_values)
		col4 = astropy.io.fits.Column(name='FLUX', format='E', array=fluxes)
		cols = astropy.io.fits.ColDefs([col1, col2, col3, col4])	
		tbhdu =astropy.io.fits.new_table(cols)
	
		prihdr = astropy.io.fits.Header()
		prihdr['COMMENT'] = "This file created by postprocessor.py from the Ultracam pipeline."
		prihdu = astropy.io.fits.PrimaryHDU(header=prihdr)
		thdulist = astropy.io.fits.HDUList([prihdu, tbhdu])
		thdulist.writeto(FITSFilename, clobber=True)



	""" For each object in the red channel try to find a match in the other two channels
	"""
	
	debug.write("Checking the red objects...", level=2)
	
	colourObjectList = []
	
	for object in redObjects:
		debug.write(str(object))
		meanPosition = object.meanPosition
		debug.write("Mean position: " + str(meanPosition))
		newIDNumber = ultracamutils.getUniqueID(colourObjectList)
		colourObject = classes.combined3ColourObject(newIDNumber)
		colourObject.setRedObject(object)

		smallestDistance = 1000
		nearestObject = greenObjects[0]
		for g in greenObjects:
			distance = ultracamutils.measureDistance(g.meanPosition, meanPosition)
			if distance < smallestDistance: 
				smallestDistance = distance
				nearestObject = g
		debug.write("Most likely match is %s at a distance of %f"%(str(nearestObject),smallestDistance), level=3)
		if (smallestDistance>float(config.MINPIXELDISTANCE)):
			debug.write("Match rejected... too far apart!")
		else: 
			colourObject.setGreenObject(nearestObject)

		smallestDistance = 1000
		nearestObject = blueObjects[0]
		for b in blueObjects:
			distance = ultracamutils.measureDistance(b.meanPosition, meanPosition)
			if distance < smallestDistance: 
				smallestDistance = distance
				nearestObject = b
		debug.write("Most likely match is %s at a distance of %f"%(str(nearestObject),smallestDistance), level=3)
		if (smallestDistance>float(config.MINPIXELDISTANCE)):
			debug.write("Match rejected... too far apart!")
		else:
			colourObject.setBlueObject(nearestObject)

		colourObjectList.append(colourObject)		
	
		
if (int(config.WRITE_JSON)==1):
	allObjects = []
	for c in colourObjectList:
		#debug.write(str(c), level = 1)
		allObjects.append(c.toJSON())

	outputFilename = ultracamutils.addPaths(config.SITE_PATH,runName) 
	outputFilename+= "_rgb.json"
	debug.write("Writing master JSON file: " + outputFilename, level = 2)

	outputfile = open( outputFilename, "w" )
	json.dump(allObjects, outputfile)
	outputfile.close()
	
	allObjects = []

	for m in redObjects:
		allObjects.append(m.toJSON())
	
	outputFilename = ultracamutils.addPaths(config.SITE_PATH,runName) 
	outputFilename+= "_r.json"
	debug.write("Writing red JSON file: " + outputFilename, level = 2)
	
	outputfile = open( outputFilename, "w" )
	json.dump(allObjects, outputfile)
	outputfile.close()


	allObjects = []
	for m in greenObjects:
		allObjects.append(m.toJSON())
	
	outputFilename = ultracamutils.addPaths(config.SITE_PATH,runName) 
	outputFilename+= "_g.json"
	debug.write("Writing green JSON file: " + outputFilename, level = 2)
	
	outputfile = open( outputFilename, "w" )
	json.dump(allObjects, outputfile)
	outputfile.close()

	allObjects = []
	for m in blueObjects:
		allObjects.append(m.toJSON())
	
	outputFilename = ultracamutils.addPaths(config.SITE_PATH,runName) 
	outputFilename+= "_b.json"
	debug.write("Writing blue JSON file: " + outputFilename, level = 2)
	
	outputfile = open( outputFilename, "w" )
	json.dump(allObjects, outputfile)
	outputfile.close()
	
"""
print "Writing a log file... just temp file..."
outfile = open("temp.log", "w")
	
MJDs = redObjects[0].getMJDs()
print MJDs
frameCounter = 0
outString = "Frame, MJD, exposure, CountsRed1, CountsRed2, CountsGreen1, CountsGreen2, CountsBlue1, CountsBlue2\n"
outfile.write(outString)
for m in MJDs:
	outString = str(frameCounter) + ", " + str(m) + ", 0.04"
	for c in colourObjectList:
		red = c.r;
		green = c.g;
		blue = c.b;
		countsRed = red.getCountsForMJD(m)
		countsGreen = green.getCountsForMJD(m)
		countsBlue = blue.getCountsForMJD(m)
		frameCounter+=1 
		outString+= ", " + str(countsRed) + ", " + str(countsGreen) + ", " + str(countsBlue)
	outString+= "\n"
	outfile.write(outString)
	
outfile.close()
"""

if (int(config.WRITE_FITS)==1):
	
	testObject = colourObjectList[8]

	redObject = testObject.r
	redObservations = redObject.getData()
	MJDs = np.array(redObservations[0])
	counts = np.array(redObservations[1])
	col1 = fits.Column(name='MJDR', format='D', array=MJDs)
	col2 = fits.Column(name='CountsR', format='D', array=counts)

	greenObject = testObject.g
	greenObservations = greenObject.getData()
	MJDs = np.array(greenObservations[0])
	counts = np.array(greenObservations[1])
	col3 = fits.Column(name='MJDG', format='D', array=MJDs)
	col4 = fits.Column(name='CountsG', format='D', array=counts)

	blueObject = testObject.b
	blueObservations = blueObject.getData()
	MJDs = np.array(blueObservations[0])
	counts = np.array(blueObservations[1])
	col5 = fits.Column(name='MJDB', format='D', array=MJDs)
	col6 = fits.Column(name='CountsB', format='D', array=counts)


	cols = fits.ColDefs([col1, col2, col3, col4, col5, col6])
	tbhdu = fits.new_table(cols)
	prihdr = fits.Header()
	prihdr['COMMENT'] = "Here's some commentary about this FITS file."
	prihdu = fits.PrimaryHDU(header=prihdr)
	thdulist = fits.HDUList([prihdu, tbhdu])
	thdulist.writeto('table.fits', clobber=True)
