#!/usr/bin/env python
import rashley_utils as utils
import sys, subprocess, re, json
import classes, numpy as np
from astropy.io import fits

	
if (__name__ == "__main__"):
	
	config = utils.readConfigFile()
	debug = classes.debugObject(config.DEBUG)

	if (len(sys.argv) < 2):
		print "Please give me a run name."
		sys.exit()
	
	runName = sys.argv[1]
	
	redFilename = utils.addPaths(config.SITE_PATH, runName) + "r.json"	
	redObjects = utils.buildObjectsFromJSON(redFilename)
	
	greenFilename = utils.addPaths(config.SITE_PATH, runName) + "g.json"
	greenObjects = utils.buildObjectsFromJSON(greenFilename)
	
	blueFilename = utils.addPaths(config.SITE_PATH, runName) + "b.json"
	blueObjects = utils.buildObjectsFromJSON(blueFilename)
		
	redObjects = utils.filterOutCosmicRays(redObjects)
	greenObjects = utils.filterOutCosmicRays(greenObjects)
	blueObjects = utils.filterOutCosmicRays(blueObjects)

	redObjects = utils.filterOutLowFrameCountObjects(redObjects, 10)
	greenObjects = utils.filterOutLowFrameCountObjects(greenObjects, 10)
	blueObjects = utils.filterOutLowFrameCountObjects(blueObjects, 10)
	
	debug.write("%d red objects after cosmic ray and low frames filtering"%(len(redObjects)))
	debug.write("%d green objects after cosmic ray and low frames filtering"%(len(greenObjects)))
	debug.write("%d blue objects after cosmic ray and low frames filtering"%(len(blueObjects)))

	print "Timing check red"
	redObjects = utils.filterOutBadTimingFrames(redObjects)
	print "Timing check green"
	greenObjects = utils.filterOutBadTimingFrames(greenObjects)
	print "Timing check blue"
	blueObjects = utils.filterOutBadTimingFrames(blueObjects)
	
		
	for i in redObjects:
		meanPosition = i.calculateMeanPosition()
	for i in greenObjects:
		meanPosition = i.calculateMeanPosition()
	for i in blueObjects:
		meanPosition = i.calculateMeanPosition()
	
	""" For each object in the red channel try to find a match in the other two channels
	"""
	
	debug.write("Checking the red objects...", level=2)
	
	colourObjectList = []
	
	for object in redObjects:
		debug.write(str(object), level=2)
		meanPosition = object.meanPosition
		debug.write("Mean position: " + str(meanPosition))
		newIDNumber = utils.getUniqueID(colourObjectList)
		colourObject = classes.combined3ColourObject(newIDNumber)
		colourObject.setRedObject(object)

		smallestDistance = 1000
		nearestObject = greenObjects[0]
		for g in greenObjects:
			distance = utils.measureDistance(g.meanPosition, meanPosition)
			if distance < smallestDistance: 
				smallestDistance = distance
				nearestObject = g
		debug.write("Most likely match is %s at a distance of %f"%(str(nearestObject),smallestDistance), level=2)
		if (smallestDistance>float(config.MINPIXELDISTANCE)):
			debug.write("Match rejected... too far apart!")
		else: 
			colourObject.setGreenObject(nearestObject)

		smallestDistance = 1000
		nearestObject = blueObjects[0]
		for b in blueObjects:
			distance = utils.measureDistance(b.meanPosition, meanPosition)
			if distance < smallestDistance: 
				smallestDistance = distance
				nearestObject = b
		debug.write("Most likely match is %s at a distance of %f"%(str(nearestObject),smallestDistance), level=2)
		if (smallestDistance>float(config.MINPIXELDISTANCE)):
			debug.write("Match rejected... too far apart!")
		else:
			colourObject.setBlueObject(nearestObject)

		colourObjectList.append(colourObject)		
	
		
if (int(config.WRITE_JSON)==1):
	allObjects = []
	for c in colourObjectList:
		debug.write(str(c), level = 1)
		allObjects.append(c.toJSON())

	outputFilename = utils.addPaths(config.SITE_PATH,runName) 
	outputFilename+= "rgb.json"
	debug.write("Writing master JSON file: " + outputFilename)

	outputfile = open( outputFilename, "w" )
	json.dump(allObjects, outputfile)
	outputfile.close()
	
	allObjects = []

	for m in redObjects:
		allObjects.append(m.toJSON())
	
	outputFilename = utils.addPaths(config.SITE_PATH,runName) 
	outputFilename+= "r.json"
	debug.write("Writing red JSON file: " + outputFilename)
	
	outputfile = open( outputFilename, "w" )
	json.dump(allObjects, outputfile)
	outputfile.close()


	allObjects = []
	for m in greenObjects:
		allObjects.append(m.toJSON())
	
	outputFilename = utils.addPaths(config.SITE_PATH,runName) 
	outputFilename+= "g.json"
	debug.write("Writing green JSON file: " + outputFilename)
	
	outputfile = open( outputFilename, "w" )
	json.dump(allObjects, outputfile)
	outputfile.close()

	allObjects = []
	for m in blueObjects:
		allObjects.append(m.toJSON())
	
	outputFilename = utils.addPaths(config.SITE_PATH,runName) 
	outputFilename+= "b.json"
	debug.write("Writing blue JSON file: " + outputFilename)
	
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
