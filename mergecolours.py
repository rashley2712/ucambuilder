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
	
	for c in channels:
		objects = allObjects[c]
		for o in objects:
			o.calculateMeanFlux()
			
	
	redObjects = allObjects['r']
	greenObjects = allObjects['g']
	for r in redObjects:
		greenDistance = 1
		closestGreen = {}
		redCoords = (r.ra, r.dec)
		for g in greenObjects:
			greenCoords = (g.ra, g.dec)
			distance = ultracamutils.calculateDistance(redCoords, greenCoords)
			if distance < greenDistance:
				closestGreen = g
				greenDistance = distance
		
		print r, closestGreen, greenDistance
		
		
	