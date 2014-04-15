#!/usr/bin/env python
import ultracamutils
import sys, subprocess, re, json
import classes, numpy as np
import astropy.io.fits
import astropy.wcs
import argparse, os, copy
import wcsclasses
	
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
		jsonFilename = ultracamutils.addPaths(config.SITE_PATH, runName) + "_" + c + ".json"	
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
	
	""" Calculate the mean flux for each object in preparation for sorting them
	"""
	for c in channels:
		objects = allObjects[c]
		for o in objects:
			o.calculateMeanFlux()
			
	""" Sort the objects
	"""
	for c in channels:
		objects = allObjects[c]
		sortedObjects = sorted(objects, key= lambda p:p.meanFlux, reverse=True)
		allObjects[c] = sortedObjects
			
	
	redObjects = allObjects['r']
	greenObjects = allObjects['g']
	blueObjects = allObjects['b']
	distanceThresholdSeconds = 10.
	distanceThresholdDegrees = distanceThresholdSeconds / 3600.

	allObjectsCopy = copy.copy(allObjects)
	
	""" Start with the colour that has the fewest number of objects
	"""
	startCatalogIndex = 'r'
	startCatalogLength = len(allObjectsCopy[startCatalogIndex])
	for c in channels:
		length = len(allObjectsCopy[c])
		if length<startCatalogLength:
			startCatalogIndex = c
			startCatalogLength = length
	debug.write("Channel with the smallest number of objects is %s with %d objects."%(channelDescriptions[startCatalogIndex], startCatalogLength))
			
	colourObjects = []
	
	threeColours = ['r', 'g', 'b']
	currentColour = startCatalogIndex
	debug.write("Running through the %s objects."%channelDescriptions[currentColour])
	threeColours.pop(threeColours.index(currentColour))
	print "Remaining colours", threeColours
	firstColourObjects = allObjectsCopy[startCatalogIndex]
	firstColourList = copy.copy(firstColourObjects)
	for o in firstColourObjects:
		id = ultracamutils.getUniqueID(colourObjects)
		colourObject = classes.combined3ColourObject(id)
		colourObject.setColourID(currentColour, o.id)
		firstColourList.remove(o)
		allObjectsCopy[currentColour] = firstColourList
		originalCoords = (o.ra, o.dec)
		for othercolour in threeColours:
			objects = allObjectsCopy[othercolour]
			debug.write("Looking for %s objects to match the current %s object"%(channelDescriptions[othercolour], channelDescriptions[currentColour]))
			closestDistance = 1
			closestObject = None
			for p in objects:
				objectCoords = (p.ra, p.dec)
				distance = ultracamutils.calculateDistance(originalCoords, objectCoords)
				if distance < closestDistance:
					closestDistance = distance
					closestObject = p
			if closestDistance < distanceThresholdDegrees:
				colourObject.setColourID(othercolour, closestObject.id)
				objects.remove(closestObject)
				allObjectsCopy[othercolour] = objects
			else: 
				debug.write("No match... to far from threshold")
			
		debug.write(colourObject)
		colourObjects.append(colourObject)
		
	debug.write("Objects remaining:")	
	for c in channels:
		debug.write("%s: %d"%(channelDescriptions[c], len(allObjectsCopy[c])))

	startCatalogIndex = threeColours[0]
	startCatalogLength = len(allObjectsCopy[startCatalogIndex])
	for c in threeColours:
		length = len(allObjectsCopy[c])
		if length<startCatalogLength:
			startCatalogIndex = c
			startCatalogLength = length
	currentColour = startCatalogIndex
			
	debug.write("Channel with the smallest number of objects is %s with %d objects."%(channelDescriptions[startCatalogIndex], startCatalogLength))
	threeColours.pop(threeColours.index(currentColour))
	print "Remaining colours:", threeColours
	
	secondColourObjects = allObjectsCopy[currentColour]
	secondColourList = copy.copy(secondColourObjects)
	for o in secondColourObjects:
		id = ultracamutils.getUniqueID(colourObjects)
		colourObject = classes.combined3ColourObject(id)
		colourObject.setColourID(currentColour, o.id)
		secondColourList.remove(o)
		allObjectsCopy[currentColour] = secondColourList
		originalCoords = (o.ra, o.dec)
		for othercolour in threeColours:
			objects = allObjectsCopy[othercolour]
			debug.write("Looking for %s objects to match the current %s object"%(channelDescriptions[othercolour], channelDescriptions[currentColour]))
			closestDistance = 1
			closestObject = None
			for p in objects:
				objectCoords = (p.ra, p.dec)
				distance = ultracamutils.calculateDistance(originalCoords, objectCoords)
				if distance < closestDistance:
					closestDistance = distance
					closestObject = p
			if closestDistance < distanceThresholdDegrees:
				colourObject.setColourID(othercolour, closestObject.id)
				objects.remove(closestObject)
				allObjectsCopy[othercolour] = objects
			else: 
				debug.write("No match... to far from threshold")
			
		debug.write(colourObject)
		colourObjects.append(colourObject)
		
	debug.write("Objects remaining:")	
	for c in channels:
		debug.write("%s: %d"%(channelDescriptions[c], len(allObjectsCopy[c])))

	currentColour = threeColours[0]
	finalObjects = allObjectsCopy[currentColour]
	finalColourList = copy.copy(finalObjects)
	for o in finalObjects:
		id = ultracamutils.getUniqueID(colourObjects)
		colourObject = classes.combined3ColourObject(id)
		colourObject.setColourID(currentColour, o.id)
		debug.write(colourObject)
		colourObjects.append(colourObject)
		finalColourList.remove(o)
		allObjectsCopy[currentColour] = finalColourList
		
	debug.write("Objects remaining:")	
	for c in channels:
		debug.write("%s: %d"%(channelDescriptions[c], len(allObjectsCopy[c])))

	
	""" Look for duplicates
	"""
	redList = []
	greenList = []
	blueList = []
	duplicates = 0
	for c in colourObjects:
		print c
	
		redID = c.getColourID('r')
		try: 
			if redList.index(redID): duplicates+=1
		except ValueError:
			if redID!= -1: redList.append(redID)

		greenID = c.getColourID('g')
		try: 
			if greenList.index(greenID): duplicates+=1
		except ValueError:
			if greenID!= -1: greenList.append(greenID)
			
		blueID = c.getColourID('b')
		try: 
			if blueList.index(blueID): duplicates+=1
		except ValueError:
			if blueID!= -1: blueList.append(blueID)
	
	print redList, greenList, blueList

	print "Duplicates: ", duplicates
	
	for c in colourObjects:
		print c.toJSON()

