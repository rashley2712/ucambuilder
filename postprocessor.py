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
	parser.add_argument('-d', '--debuglevel', default = 1, type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	arg = parser.parse_args()
	
	config = ultracamutils.readConfigFile(arg.configfile)
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);
	debug.write(config)
	
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
		jsonFilename = ultracamutils.addPaths(config.WORKINGDIR, runName) + "_" + c + "_raw.json"	
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
			debug.write("%d %s objects (was...%d) after cosmic ray filtering"%(len(allObjects[c]), channelDescriptions[c], beforeCount), level =2 )
	
	percentThreshold = 20
	for c in channels:
			objects = allObjects[c]
			beforeCount = len(objects)
			objects = ultracamutils.filterOutLowFrameCountObjects(objects, percentThreshold)
			allObjects[c] = objects
			debug.write("%d %s objects (was...%d) after removing objects that appear on fewer than %d%% of the frames"%(len(allObjects[c]), channelDescriptions[c], beforeCount, percentThreshold), level = 2)
	 
	pixelThreshold = 1.
	for c in channels:
			objects = allObjects[c]
			beforeCount = len(objects)
			for o in objects:
				meanFWHM = o.calculateMeanFWHM()
			objects = ultracamutils.filterOutPixels(objects, pixelThreshold)
			allObjects[c] = objects
			debug.write("%d %s objects (was...%d) after removing objects that have a pixel size smaller than %d."%(len(allObjects[c]), channelDescriptions[c], beforeCount, pixelThreshold), level = 2)

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
			debug.write(i)
	
	if (arg.xyls):
		for c in channels:
			debug.write("Creating an XYLS file for the %s channel"%(channelDescriptions[c]))
			objects = allObjects[c]
			sortedObjects = sorted(objects, key= lambda p:p.meanFlux, reverse=True)
			x_values, IDs, y_values, fluxes = [], [], [], []
		
			maxFlux = sortedObjects[1].meanFlux
			fluxThreshold = maxFlux/10000.
			debug.write("Max flux:" + str(maxFlux) + "Theshold:" + str(fluxThreshold))
			
			for num, i in enumerate(sortedObjects):
				if i.meanFlux<fluxThreshold: break;
				IDs.append(i.id)
				x_values.append(i.meanPosition[0])
				y_values.append(i.meanPosition[1])
				fluxes.append(i.meanFlux)

			runDate, runNumber = ultracamutils.separateRunNameAndDate(arg.runname)
			FITSFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_" + c + ".xyls"
			debug.write("Writing FITS file: " + FITSFilename, level=2)
	
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
	"""
		
if (int(config.WRITE_JSON)==1):
	for c in channels:
		outputFilename = ultracamutils.addPaths(config.SITE_PATH,runName) 
		outputFilename+= "_" + c + ".json"
		debug.write("Writing the refined %s object catalogs to .json file: %s"%(channelDescriptions[c], outputFilename), level =2 )
		objects = allObjects[c]
		sortedObjects = sorted(objects, key= lambda p:p.meanFlux, reverse=True)
		
		jsonArray = []
		for m in objects:
			jsonArray.append(m.toJSON())
	
	
		outputfile = open( outputFilename, "w" )
		json.dump(jsonArray, outputfile)
		outputfile.close()


