#!/usr/bin/env python
import ultracamutils, ucamObjectClass
import sys, subprocess, re, json
import classes, numpy as np
import astropy.io.fits
import astropy.wcs
import argparse, os, copy
import wcsclasses
	
	
def addPhotometry(colourObject, colour, exposureArray):
	exposureArray = o.exposures
	for e in exposureArray:
		# Match the exposure with a frame based on the MJD
		MJD = e.MJD
		frameIndex = -1
		for frame in frameData:
			if frame.MJD == MJD:
				frameIndex = frame.frameIndex
				frameObject = frame
				break
			
		newExposure = {'frameIndex': frameIndex}
		newExposure['magnitude'] = e.counts
		newExposure['fwhm'] = e.FWHM
		newExposure['position'] = e.centroid
			
		colourObject.addExposure(colour, newExposure)
		

	

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
	pixelMatch = False      # This is set to True if we don't have a WCS solution and need to use pixel locations for matching
	
	""" Load the information about the frames
	"""
	jsonFilename = ultracamutils.addPaths(config.SITE_PATH, runName) + "_frameInfo.json"
	debug.write("Loading frame info from %s"%(jsonFilename), level = 2)
	jsonFile = open(jsonFilename, 'r')
	jsonObjects = json.loads(jsonFile.read())
	
	frameData = []
	
	for j in jsonObjects:
		object = json.loads(j)
		frame = ucamObjectClass.frameObject()
		frame.setFromObject(object)
		frameData.append(frame)

	debug.write("Loaded info for %d frames."%(len(frameData)))
	
	""" Load the objects from the .json files.... channel by channel (r, g, b)
	"""
	for c in channels:
		jsonFilename = ultracamutils.addPaths(config.WORKINGDIR, runName) + "_" + c + ".json"	
		debug.write("Loading the json file for the %s objects from path: %s"%(channelDescriptions[c], jsonFilename), level = 2)
		objects = ultracamutils.buildObjectsFromJSON(jsonFilename)
		allObjects[c] = objects
		debug.write("%d %s objects loaded."%(len(allObjects[c]), channelDescriptions[c]), level = 2)
		
	""" Look for and load WCS solutions for each colour (if they exist)
	"""
	for c in channels:
		wcsSolutionFilename = ultracamutils.addPaths(config.WORKINGDIR, runName) + "_" + c + ".wcs"
		if os.path.exists(wcsSolutionFilename):
			debug.write("There is a WCS solution for channel: %s"%(channelDescriptions[c]), level = 2)
			wcsParametersFile = astropy.io.fits.open(wcsSolutionFilename)
			header = wcsParametersFile[0].header
			wcs = wcsclasses.wcsSolution()
			
			equinox = float(header['EQUINOX'])
			referenceCoord = float(header['CRVAL1']), float(header['CRVAL2'])
			referencePixel = float(header['CRPIX1']), float(header['CRPIX2'])
			
			CD_array = [ [header['CD1_1'], header['CD1_2']], [ header['CD2_1'], header['CD2_2'] ] ]
			
			wcs.setSolution(equinox, referenceCoord, referencePixel, CD_array)

			w = astropy.wcs.WCS(wcsParametersFile[0].header)

			wcsParametersFile.close()

			for o in allObjects[c]:
				(x, y) = o.calculateMeanPosition()
				(ra, dec) = w.all_pix2world(x,y, 1.)
				o.setWorldPosition(ra, dec)			
		else:
			debug.write("No WCS solution... will have to fall back to 'pixel' matching.", level = 2)
			pixelMatch = True
	
	""" Calculate the mean flux for each object in preparation for sorting them. Also calculate their mean pixel position.
	"""
	for c in channels:
		objects = allObjects[c]
		for o in objects:
			o.calculateMeanFlux()
			o.calculateMeanPosition()
			
	""" Sort the objects
	"""
	for c in channels:
		objects = allObjects[c]
		sortedObjects = sorted(objects, key= lambda p:p.meanFlux, reverse=True)
		allObjects[c] = sortedObjects
			
	
			
	colour = 'r'
	objects = allObjects[colour]
	masterObjectList = []
	for o in objects:
		newIDNumber = ultracamutils.getUniqueID(masterObjectList)
		colourObject = ucamObjectClass.colourObject(newIDNumber)
		debug.write("Created a new colourObject with id: %d"%(newIDNumber))
		
		colourObject.setMeanPosition(colour, o.meanPosition)
		colourObject.colourID[colour] = o.id
		
		""" Now move the photometry into the new object
		"""
		addPhotometry(colourObject, colour, o.exposures)
		masterObjectList.append(colourObject)

	distanceThreshold = float(config.MINPIXELDISTANCE)
	print "Threshold", distanceThreshold

	colour = 'g'
	objects = allObjects[colour]
	for o in objects:
		""" First see if we have a position match in our existing objects
		"""
		position = o.meanPosition
		closestDistance = 1000
		closestObject = None
		for m in masterObjectList:
			r_distance = ultracamutils.calculateDistance(position, m.meanPosition['r'])
			if r_distance < closestDistance: 
				closestDistance = r_distance
				closestObject = m
		if closestDistance > distanceThreshold:
			closestObject = None   # Reject the match if it is too far away
		
		if closestObject == None:
			newIDNumber = ultracamutils.getUniqueID(masterObjectList)
			colourObject = ucamObjectClass.colourObject(newIDNumber)
			debug.write("Could find no match to this green object!")
			debug.write("Created a new colourObject with id: %d"%(newIDNumber))
			colourObject.setMeanPosition(colour, o.meanPosition)
			colourObject.colourID[colour] = o.id
			addPhotometry(colourObject, colour, o.exposures)
			masterObjectList.append(colourObject)
		else: 
			closestObject.setMeanPosition(colour, o.meanPosition)
			closestObject.colourID[colour] = o.id
		
			addPhotometry(closestObject, colour, o.exposures)
			
	colour = 'b'
	objects = allObjects[colour]
	for o in objects:
		""" First see if we have a position match in our existing objects
		"""
		position = o.meanPosition
		closestDistance = 1000
		closestObject = None
		for m in masterObjectList:
			if m.colourID['r']!= -1 :
				r_distance = ultracamutils.calculateDistance(position, m.meanPosition['r'])
				if r_distance < closestDistance: 
					closestDistance = r_distance
					closestObject = m
			if m.colourID['g']!= -1 :
				g_distance = ultracamutils.calculateDistance(position, m.meanPosition['g'])
				if g_distance < closestDistance:
					closestDistance = g_distance
					closestObject = m
					
		if closestDistance > distanceThreshold:
			closestObject = None   # Reject the match if it is too far away
		
		
		if closestObject == None:
			newIDNumber = ultracamutils.getUniqueID(masterObjectList)
			colourObject = ucamObjectClass.colourObject(newIDNumber)
			debug.write("Could find no match to this blue object!")
			debug.write("Created a new colourObject with id: %d"%(newIDNumber))
			colourObject.setMeanPosition(colour, o.meanPosition)
			addPhotometry(colourObject, colour, o.exposures)
			colourObject.colourID[colour] = o.id
			masterObjectList.append(colourObject)
		else: 
			closestObject.setMeanPosition(colour, o.meanPosition)
			closestObject.colourID[colour] = o.id
			addPhotometry(closestObject, colour, o.exposures)
			
		
	""" Now write out the objectInfo to a JSON file...
	"""
	outputFilename = ultracamutils.addPaths(config.SITE_PATH, arg.runname)
	outputFilename+= "_objects.json"

	debug.write("Writing %d objects to: %s"%(len(masterObjectList), outputFilename))
		
	JSONObjects = []
	
	for m in masterObjectList:
		JSONObjects.append(m.toJSON())
	
	outputfile = open(outputFilename, "w")
	
	json.dump(JSONObjects, outputfile)
	outputfile.close()

	
