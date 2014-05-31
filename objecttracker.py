#!/usr/bin/env python
import astropy, sys, time, os, subprocess, tempfile
import astropy.io.fits
from scipy import *
import numpy, json
from trm import ultracam
from trm.ultracam.UErrors import PowerOnOffError, UendError, UltracamError
import rashley_utils as utils
import classes
import Image,ImageDraw

def getNextFrame():
	global rdat, frameCounter, actualFrame, frameMJD, nxmax, nymax, numWindows, frameInfo, goodTime
	try:
		ccdFrame = rdat()
	except UltracamError:
		print "There was an error reading the next frame."
		print UltracamError
		return

	frameR = ccdFrame[0]
	frameG = ccdFrame[1]
	frameB = ccdFrame[2]

	if (frameCounter == 1):
		# Work out the window structure for this run...
		nxmax = frameR.nxmax
		nymax = frameR.nymax
		numWindows = len(frameR)
		debug.write("Frame dimensions: (%d, %d), Num windows %d"%(nxmax, nymax, numWindows), level = 2)
		frameInfo.channel = channel
		frameInfo.nxmax = nxmax
		frameInfo.nymax = nymax
		for nwin,win in enumerate(frameR._data):
			frameInfo.addWindow(win.llx, win.lly, win.nx, win.ny)
		debug.write(frameInfo, level = 2)

   	frameMJD = frameR.time.mjd
	goodTime = frameR.time.good 
   	#timeTest = rdat.time()
	#timeTest = rtime()
	#debug.write(timeTest)
	debug.write("Frame: " + str(frameCounter) + "[" + str(actualFrame) + "]" + "  MJD: " + str(frameMJD), level = 1)
	
	frameWindows = []
	if (channel == 'r'):
		for i in frameR._data:
			frameWindows.append(i.data)
	if (channel == 'g'):
		for i in frameG._data:
			frameWindows.append(i.data)
	if (channel == 'b'):
		for i in frameB._data:
			frameWindows.append(i.data)
		
	return frameWindows
	
def updateMasterObjectList(objects, windowObject):
	global frameMJD, masterObjectList
	
	for ob in objects:
		# add x, y offset for this window
		debug.write("sextractor x, y (%d, %d)"%(int(ob['x']),int(ob['y'])), level = 3)
		ob['x'] = ob['x'] + windowObject.xll
		ob['y'] = ob['y'] + windowObject.yll
		debug.write("actual x, y (%d, %d)"%(int(ob['x']),int(ob['y'])), level = 3)
		# is this object in the master list yet?
		objectRecognised = False
		for eo in masterObjectList:
			if (eo.isDistanceMatch(ob)!=-1):
				eo.addExposureByObject(ob, frameMJD)
				objectRecognised = True
				break
		if (not objectRecognised):
			newIDNumber = utils.getUniqueID(masterObjectList)
			debug.write("New object detected... assigning a new ID number: %d"%newIDNumber)
			newObject = classes.ObservedObject(newIDNumber)
			newObject.addExposureByObject(ob, frameMJD)
			masterObjectList.append(newObject)

config = utils.readConfigFile()
debug = classes.debugObject(config.DEBUG)

if (len(sys.argv) < 2):
	print "Please give me a run name."
	sys.exit()

runName = sys.argv[1]
runFilename = utils.addPaths(config.ULTRACAMRAW, runName)

channel = 'r'

startFrame = 1
requestedNumFrames = 10
keepTmpFiles = False
frameInfo = classes.FrameObject()

for i in sys.argv[2:]:
	if i[:2]=="-n": 
		print "number of frames:",  i[2:] 
		framesStr = i[2:]
		requestedNumFrames = int(framesStr)
	if i[:2]=="-s": 
		print "starting at frame:", i[2:]
		framesStr = i[2:]
		startFrame = int(framesStr)
	if i[:2]=="-o": 
		print "output file:",  i[2:] 
		outputFilename = i[2:]
	if i[:2]=="-k":
		keepTmpFiles = True
	if i[:2]=="-c":
		channel = i[2:]
	if i[:2]=="-d":
		debugLevel = i[2:]
		debug.setLevel(debugLevel)

outputFilename = utils.addPaths(config.SITE_PATH,runName)  + "deepimage" + channel + "."
debug.toggleTimeLog()

rdat  = ultracam.Rdata(runFilename, startFrame, server=False)
numFrames = rdat.ntotal()
if (requestedNumFrames==-1): requestedNumFrames = numFrames;
frameMJD = 0

debug.write("Total frames in %s: %d"%(runName,numFrames))
debug.write("nBlue: " + str(rdat.nblue))

masterObjectList = []
summedFrames = []
frameCounter = 0
firstFrame = True
goodTime = False
for i in range(requestedNumFrames):
	actualFrame = frameCounter + startFrame 
	if (actualFrame<numFrames):
		frameCounter+= 1
	else: 
		break
		
	frameWindows = getNextFrame()

	if frameWindows == None:
		print "We couldn't read the next frame... exiting the loop"
		break
		
	if goodTime == False:
		print "Timing error.... skipping the frame"
		continue

	if (channel == 'b') & (frameCounter % rdat.nblue != 0):      # This is an empty blue frame so skip it.
		continue   
	
	if (firstFrame):	
		for j in range(frameInfo.numWindows):
			width = frameInfo.getWindow(j).xsize
			height = frameInfo.getWindow(j).ysize
			summedFrame = numpy.zeros((width, height))
			summedFrames.append(summedFrame)
		firstFrame = False
			
	for j in range(frameInfo.numWindows): 
	
		frameImage = frameWindows[j]
		utils.createFITS(frameCounter, frameImage)
		utils.runSex(frameCounter)
		
		frameImage = numpy.rot90(frameImage)
							
		summedFrame = summedFrames[j]
		summedFrame = numpy.add(summedFrame, frameImage)
		summedFrames[j] = summedFrame
			
		if (int(config.KEEP_TMP_FILES)==0):utils.removeFITS(frameCounter);
		
		objects = utils.readSexObjects(frameCounter)
		
		if (int(config.KEEP_TMP_FILES)==0):utils.removeCAT(frameCounter)
		
		objects = utils.rejectBadObjects(objects)
		
		updateMasterObjectList(objects, frameInfo.getWindow(j))
	
		debug.write("Total objects in this window:" + str(len(objects)))


	
averageFrames = []	
normalisedImages = []
for j in range(frameInfo.numWindows):
	averageFrame = numpy.divide(summedFrames[j], frameCounter)
	averageFrames.append(averageFrame)
	normalisedImage = utils.percentiles(averageFrame, 25, 99)
	normalisedImages.append(normalisedImage)

	
# Construct full frame from Windows
fullFrame = numpy.zeros((frameInfo.nxmax, frameInfo.nymax))
for j in range(frameInfo.numWindows):
	xll = frameInfo.getWindow(j).xll 
	yll = frameInfo.getWindow(j).yll 
	xsize = frameInfo.getWindow(j).xsize 
	ysize = frameInfo.getWindow(j).ysize 
	rotatedImage = numpy.flipud(normalisedImages[j])
	fullFrame[xll:xll+xsize, yll:yll+ysize] = fullFrame[xll:xll+xsize, yll:yll+ysize] + rotatedImage
	
fullFrame = numpy.fliplr(fullFrame)
	
debug.write("Objects detected:")
for m in masterObjectList:
	debug.write(m.toString())


imgData = fullFrame
imgSize = numpy.shape(imgData)
imgLength = imgSize[0] * imgSize[1]
testData = numpy.reshape(imgData, imgLength, order="F")
img = Image.new("L", imgSize)
img.putdata(testData)
debug.write("Writing PNG file: " + outputFilename + ".png") 
img.save(outputFilename + "png", "PNG")


if (int(config.WRITE_FITS)==1): 
	FITSFilename = outputFilename + "fits"
	debug.write("Writing FITS file: " + FITSFilename, level=1)
	utils.saveFITSImage(fullFrame, FITSFilename)


if (int(config.WRITE_JSON)==1):

	allObjects = []

	for m in masterObjectList:
		allObjects.append(m.toJSON())
		print m
	
	outputFilename = utils.addPaths(config.SITE_PATH,runName) 
	outputFilename+= channel + ".json"
	debug.write("Writing JSON file: " + outputFilename)
	
	outputfile = open( outputFilename, "w" )
	json.dump(allObjects, outputfile)
	outputfile.close()
	
	
