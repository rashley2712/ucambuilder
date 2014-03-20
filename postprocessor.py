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
	
	runName = arg.runname

	debug = classes.debugObject(arg.debuglevel)
	debug.write("Getting run info from the file:" + config.RUNINFO, level = 2)
	runInfo = ultracamutils.getRunInfo(config.RUNINFO, arg.runname)
	
	debug.write("Run Info:\n----------------------", level = 2)
	debug.write(runInfo, level = 2)
	debug.write("----------------------", level = 2)

	channels = ['r', 'g', 'b']
	channelDescriptions = {'r': "Red", 'g': "Green", 'b': "Blue"}
	allObjects = {'r': [], 'g':[], 'b':[]}
	
	""" Load the objects from the .json files.... channel by channel (r, g, b)
	"""
	for c in channels:
		debug.write("Loading the json file for the %s objects."%(channelDescriptions[c]), level = 2)
		jsonFilename = ultracamutils.addPaths(config.SITE_PATH, runName) + "_" + c + "_raw.json"	
		objects = ultracamutils.buildObjectsFromJSON(jsonFilename)
		allObjects[c] = objects
		debug.write("%d %s objects loaded."%(len(allObjects[c]), channelDescriptions[c]), level = 2)
	
	
	""" Do some filtering of the objects
	"""
	for c in channels:
			objects = allObjects[c]
			beforeCount = len(objects)
			objects = ultracamutils.filterOutCosmicRays(objects)
			allObjects[c] = objects
			debug.write("%d %s objects (was...%d) after cosmic ray filtering"%(len(allObjects[c]), channelDescriptions[c], beforeCount))
	
	percentThreshold = 20
	for c in channels:
			objects = allObjects[c]
			beforeCount = len(objects)
			objects = ultracamutils.filterOutLowFrameCountObjects(objects, percentThreshold)
			allObjects[c] = objects
			debug.write("%d %s objects (was...%d) after removing objects that appear on fewer than %d%% of the frames"%(len(allObjects[c]), channelDescriptions[c], beforeCount, percentThreshold))
	 
	pixelThreshold = 1.
	for c in channels:
			objects = allObjects[c]
			beforeCount = len(objects)
			for o in objects:
				meanFWHM = o.calculateMeanFWHM()
			objects = ultracamutils.filterOutPixels(objects, pixelThreshold)
			allObjects[c] = objects
			debug.write("%d %s objects (was...%d) after removing objects that have a pixel size smaller than %d."%(len(allObjects[c]), channelDescriptions[c], beforeCount, pixelThreshold))

	""" Perform the timing checks
		
	for c in channels:
		debug.write("Timing check %s"%(channelDescriptions[c]))
		objects = allObjects[c]
		objects = ultracamutils.filterOutBadTimingFrames(objects)
		allObjects[c] = objects
	"""
	
	
	""" Perform some calculations on the objects themselves to get meanFlux and meanPosition for each one
	"""
	for c in channels:
		objects = allObjects[c]
		for i in objects:
			meanPosition = i.calculateMeanPosition()
			meanFlux = i.calculateMeanFlux()
			meanFWHM = i.calculateMeanFWHM()
			print i
	
	if (arg.xyls):
		for c in channels:
			debug.write("Creating an XYLS file for the %s channel"%(channelDescriptions[c]))
			objects = allObjects[c]
			sortedObjects = sorted(objects, key= lambda p:p.meanFlux, reverse=True)
			x_values, IDs, y_values, fluxes = [], [], [], []
		
			maxFlux = sortedObjects[1].meanFlux
			fluxThreshold = maxFlux/1000.
			print "Max flux:", maxFlux, "Theshold:", fluxThreshold
			
			for num, i in enumerate(sortedObjects):
				if i.meanFlux<fluxThreshold: break;
				IDs.append(i.id)
				x_values.append(i.meanPosition[0])
				y_values.append(i.meanPosition[1])
				fluxes.append(i.meanFlux)

			print fluxes

			runDate, runNumber = ultracamutils.separateRunNameAndDate(arg.runname)
			FITSFilename = ultracamutils.addPaths(config.SITE_PATH, arg.runname) + "_" + c + ".xyls"
			debug.write("Writing FITS file: " + FITSFilename, level=1)
	
			col1 = astropy.io.fits.Column(name='ID', format='I', array=IDs)
			col2 = astropy.io.fits.Column(name='X', format='E', array=x_values)
			col3 = astropy.io.fits.Column(name='Y', format='E', array=y_values)
			col4 = astropy.io.fits.Column(name='FLUX', format='E', array=fluxes)
			cols = astropy.io.fits.ColDefs([col1, col2, col3, col4])	
			tbhdu =astropy.io.fits.new_table(cols)
	
			prihdr = astropy.io.fits.Header()
			prihdr['TARGET'] = runInfo.target
			prihdr['RA'] = runInfo.ra * 15.
			prihdr['DEC'] = runInfo.dec
			prihdr['COMMENT'] = "This file created by postprocessor.py from the Ultracam pipeline."
			prihdr['RUNIDENT'] = arg.runname
			prihdr['CHANNEL'] = channelDescriptions[c]
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
