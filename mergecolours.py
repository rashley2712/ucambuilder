#!/usr/bin/env python
import ultracamutils
import sys, subprocess, re, json
import classes, numpy as np
import astropy.io.fits
import astropy.wcs
import argparse, os
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
	distanceThresholdSeconds = 5.
	distanceThresholdDegrees = distanceThresholdSeconds / 3600.
	colourObjects = []
	for r in redObjects:
		id = ultracamutils.getUniqueID(colourObjects)
		colourObject = classes.combined3ColourObject(id)
		colourObject.setRedObject(r)
		greenDistance = 1
		closestGreen = {}
		redCoords = (r.ra, r.dec)
		for g in greenObjects:
			greenCoords = (g.ra, g.dec)
			distance = ultracamutils.calculateDistance(redCoords, greenCoords)
			if distance < greenDistance:
				closestGreen = g
				greenDistance = distance
				
		if (greenDistance < distanceThresholdDegrees):
			colourObject.setGreenObject(closestGreen)
			colourObject.rgDistance = greenDistance
		else: 
			print "too far for a match"
			colourObject.setGreenObject(None)
			colourObject.rgDistance = 1000

		blueDistance = 1
		for b in blueObjects:
			blueCoords = (b.ra, b.dec)
			distance = ultracamutils.calculateDistance(redCoords, blueCoords)
			if distance < blueDistance:
				closestBlue = b
				blueDistance = distance

		if (blueDistance < distanceThresholdDegrees):
			colourObject.setBlueObject(closestBlue)
			colourObject.rbDistance = blueDistance
		else: 
			print "too far for a match"
			colourObject.setBlueObject(None)
			colourObject.rbDistance = 1000
		
		greenCoords = (closestGreen.ra, closestGreen.dec)
		blueCoords = (closestBlue.ra, closestBlue.dec)
		gbDistance = ultracamutils.calculateDistance(greenCoords, blueCoords)
		colourObject.gbDistance = gbDistance

		colourObjects.append(colourObject)
		print colourObject.summaryString()
		
	""" Do some diagnostics on the top matches
	"""
	for n, c in enumerate(colourObjects):
		print "ID: ", n
		print c
		r = c.r
		g = c.g
		b = c.b
		if (r!=None): print "Red   RA:%10.5f DEC:%10.5f"%(r.ra, r.dec)
		if (g!=None): print "Green RA:%10.5f DEC:%10.5f"%(g.ra, g.dec)
		if (b!=None): print "Blue  RA:%10.5f DEC:%10.5f"%(b.ra, b.dec)
		print c.rgDistance, c.rbDistance, c.gbDistance
		if n>100: break
		
	""" Look for duplicates
	"""
	rIDlist = []
	gIDlist = []
	bIDlist = []
	for c in colourObjects:
		if (c.r!=None): 
			try:
				index = rIDlist.index(c.r.id)
				print "Red duplicate:", c.r.id, index
			except ValueError:
				rIDlist.append(c.r.id)
				print "Adding index: ", c.r.id, len(rIDlist)

		if (c.g!=None): 
			try:
				index = gIDlist.index(c.g.id)
				print "Green duplicate:", index
			except ValueError:
				gIDlist.append(c.g.id)
		if (c.b!=None): 
			try:
				index = bIDlist.index(c.b.id)
				print "Blue duplicate:", index
			except ValueError:
				bIDlist.append(c.b.id)
		
	print rIDlist
		
		
		
	
