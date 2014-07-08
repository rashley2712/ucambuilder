#!/usr/bin/env python
import ultracamutils, ucamObjectClass
import sys, subprocess, re, json, itertools
import classes, numpy, matplotlib.pyplot 
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
		
def filter3ColourObjects(objects):
	""" Returns a smaller list containing only those objects that have photometry in all three channels
	"""
	filteredList = []
	for o in objects:
		if (o.colourID['r']!=-1) & (o.colourID['g']!=-1) & (o.colourID['b']!=-1): 
			filteredList.append(o)
	
	return filteredList
	
def computeFrameOverlap(photometry1, photometry2):
	""" Compute the number of frames that overlap in 2 frame lists
	"""	
	overlapCount = 0
	for index1, measurement1 in photometry1: 
		for index2, measurement2 in photometry2:
			if index2==index1: overlapCount+= 1
		
	return overlapCount
	
def computeDifferentialPhotometry(photometry1, photometry2):
	""" Takes two arrays of data, divides one by the other... only returns values where the data arrays overlap
	"""	
	differentialPhotometry = []
	
	
	for index1, measurement1 in photometry1: 
		for index2, measurement2 in photometry2:
			if index2==index1:
				diffPhot = {'frameIndex': index1, 'diffMagnitude': measurement1/measurement2}
				differentialPhotometry.append(diffPhot)
				
	return differentialPhotometry
	
def findMaxFramesByColour(colour, objectList):
	numFrames = 0
	for o in objectList:
		photometry = o.getAllPhotometryByColour(colour)
		if len(photometry)>numFrames: numFrames = len(photometry)
	
	return numFrames
	
def getComparisonPhotometry(colour, frameIndex):
	for f in frameData:
		if f.frameIndex==frameIndex:
			return f.comparisonPhotometry[colour]
	return -1
	


if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(description='Reads the files produced by earlier steps in the pipeline.')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	parser.add_argument('--xyls', action='store_true', help='Create an XY-list for each channel. Used by Astrometry.net')
	parser.add_argument('-d', '--debuglevel', default = 1, type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('-s', '--saveasdiff', action='store_true', help='Compute differential photometry and save all photometry as fraction of comparison stars.')
	
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
	colours = ['r', 'g', 'b']
	
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

	""" Now everything is in the masterObjectList, we can delete the loaded objects """
	del allObjects


	""" Now look for comparison objects
	    Do this independently for each colour
		1. Find objects that have photometry for most frames
	"""
	for colour in channels:
		print "Switching to " + channelDescriptions[colour]
		numFrames = findMaxFramesByColour(colour, masterObjectList)
		print "max. frames for %s is %d."%(channelDescriptions[colour], numFrames)
		colourComparisons = []
		for o in masterObjectList:
			numExposures = len(o.photometry[colour])
			coverage = float(numExposures) / float(numFrames) * 100.
			print "%d has %d frames or %f%%."%(o.id, numExposures, coverage)
			potentialComparison = {'id': o.id, 'coverage': coverage}
			colourComparisons.append(potentialComparison)
		colourComparisons = sorted(colourComparisons, key= lambda c:c['coverage'], reverse=True)

		coverageThreshold = float(config.COMPARISON_THRESHOLD)
		filteredComparisons = []
		for comparison in colourComparisons:
			if comparison['coverage']>coverageThreshold:
				filteredComparisons.append(comparison)

		debug.write("Of the %d objects, %d have greater than %4.2f%% coverage in the %s channel."%(len(masterObjectList), len(filteredComparisons), coverageThreshold, channelDescriptions[colour]))

		if (len(filteredComparisons))>10:
			filteredComparisons = filteredComparisons[:10]
			debug.write("Trimmed comparisons down to %d", len(filteredComparisons))

		if len(filteredComparisons)<3:
			debug.write("Fewer than 3 objects can be used for comparisons... abandoning automatic comparison detection.", level = 2)
		else: 
			""" Compare the scatter on each....
			"""
			combinations = itertools.combinations(filteredComparisons, 2)

			compareList = []

			for i, c in enumerate(combinations):
				print "Calculating differential photometry for object: %d vs object %d"%(c[0]['id'], c[1]['id'])
				object1 = ultracamutils.getObjectByID(masterObjectList, c[0]['id'])
				object2 = ultracamutils.getObjectByID(masterObjectList, c[1]['id'])
				object1Photometry = object1.getAllPhotometryByColour('r')
				object2Photometry = object2.getAllPhotometryByColour('r')
				frameOverlap = computeFrameOverlap(object1Photometry, object2Photometry)
				overlapPercent = float(frameOverlap)/float(numFrames) * 100. 
				if (overlapPercent>coverageThreshold):
					differentialPhotometry = computeDifferentialPhotometry(object1Photometry, object2Photometry)
					print "Frame overlap:", frameOverlap, overlapPercent

					frames = [ p['frameIndex'] for p in differentialPhotometry ]
					mags = [ p['diffMagnitude'] for p in differentialPhotometry]
					mean = numpy.mean(mags)
					stddev = numpy.std(mags)
					variation = stddev / mean
					comparison = {'p1id': c[0]['id'], 'p2id': c[1]['id'], 'cov': variation}
					compareList.append(comparison)

					print "Mean", mean, "Std dev", stddev, "Co-efficient of variation", variation
				else:
					debug.write("These two objects (id1: %d and id2: %d) had less than %4.2f overlapping data points... cannot be checked... abandoning them"%(c[0]['id'], c[1]['id'], coverageThreshold), level = 2) 


				#matplotlib.pyplot.subplot()
				#matplotlib.pyplot.scatter(frames, mags)
				#matplotlib.pyplot.show()


			covs = [ c['cov'] for c in compareList]
			

			covs_mean = numpy.mean(covs)
			covs_stddev = numpy.std(covs)

			print "Mean COV", covs_mean, "STDDEV COV", covs_stddev, "len", len(covs)

			""" Sort the compareList """
			compareList = sorted(compareList, key= lambda c:c['cov'], reverse=False)

			""" Choose the pair with the lowest 'common co-efficient of variance' """
			chosenPair = compareList[0]

			print "Chosen:", chosenPair
			comparison1 = ultracamutils.getObjectByID(masterObjectList, chosenPair['p1id'])
			comparison2 = ultracamutils.getObjectByID(masterObjectList, chosenPair['p2id'])

			comparison1.setComparisonFlag(colour)
			comparison2.setComparisonFlag(colour)

			print "1:", comparison1
			print "2:", comparison2
		
			print "Computing frame comparison photometry for colour %s"%(channelDescriptions[colour])
			# For each frame in the frame list, create a comparison magnitude, equal to the average of the comparisons on that frame
			for f in frameData:
				frameIndex = f.frameIndex
				photometry1 = comparison1.getPhotometry(colour, frameIndex)
			
				if (photometry1!=-1):
					comparison = photometry1
					f.setComparisonPhotometry(colour, comparison)
				else:
					f.setComparisonPhotometry(colour, -1)
				
			
	
	#sys.exit()
			
	""" Now look for comparison objects
	    Choose the brightest object with info in all colours
	
	filteredObjects = filter3ColourObjects(masterObjectList)
	print "%d objects in total. %d have photometry in all three colours."%(len(masterObjectList), len(filteredObjects))
	
	for o in filteredObjects:
		o.calculateSigma()
		
	for c in colours:
		print c
		deviations = []
		for o in filteredObjects:
			deviations.append(o.deviations[c])
		print deviations
		stddev = numpy.std(deviations)
		mean = numpy.mean(deviations)
		print stddev
		#Which ones lie withing 1 sigma of this stddev?
		for o in filteredObjects:
			if abs(o.deviations[c] - mean) < stddev:
				print "Comparison: " + str(o)
				o.setComparisonFlag(c)
			else:
				print "Variable: " + str(o)
	
	#Check if the object passes the 1 sigma test in all three channels
	for o in filteredObjects:
		print o.comparisonFlags
		o.testComparison()
		
	for o in filteredObjects:
		if o.isComparison:
			print "Comparison: " + str(o)
		else:
			print "Variable: " + str(o)
			
	#For each frame in the frame list, create a comparison magnitude, equal to the average of the comparisons on that frame
	comparisonMagnitudes = []
	for f in frameData:
		frameIndex = f.frameIndex
		for o in filteredObjects:
			for c in colours:
				magnitude = o.getPhotometry(c, frameIndex)
				
		
	"""
	
	if arg.saveasdiff:
		for o in masterObjectList:
			print "Computing differential photometry for object %d"%(o.id)
			for c in channels:
				print channelDescriptions[c]
				if o.comparisonFlags[c]==True:
					print "skipping .... is a comparison"
				else:
					photometry = o.getAllPhotometryByColour(c)
					
					for p in photometry:
						frameIndex, magnitude = p
						comparison = getComparisonPhotometry(c, frameIndex)
						if comparison==-1:
							o.removePhotometry(c, frameIndex)
						else:
							relativePhotometry = magnitude/comparison
							o.setPhotometry(c, frameIndex, relativePhotometry)
							
					
						
					#for f in frameData:
					#	comparison = f.comparisonPhotometry[c]
					#	print comparison
				
				
	
		
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
	
	""" Also re-write the frame data
	"""
	outputFilename = ultracamutils.addPaths(config.SITE_PATH, arg.runname)
	outputFilename+= "_frameInfo.json"
		
	frameObjects = []
		
	for f in frameData:
		frameObjects.append(f.toJSON())
		
	outputfile = open(outputFilename, "w")
		
	json.dump(frameObjects, outputfile)

	outputfile.close()
		

	
